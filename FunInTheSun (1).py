# FunInTheSun.py Open Beta


#Matthew Barzal
#Nasa GSFC 2023 summer intern, University of Maryland, College Park


#This program runs /home/jovyan/efs/jrclark/PFe-Counting.ipynb as a straight .py 


#Credit to Julia Clark, USBC
#Much of the meat and potatoes of this program was orginally written by Julia Clark as a .ipynb. It was adapted for automated usage by any user as a .py program by Matthew Barzal. A number of lines were added while some were altered to prevent progress loss and enhance fluidity.


'''
read me.txt
Program description

General operating 
    This program generates a database of Polar Faculae counts from SDO's HMI instrument. It outputs monthly csv files with counts, a smoothed counts csv, and plots of the counts over time for each month.

    It is set up to only read one month at a time. In other words, when it prompts for a start date and end date, the user should only enter dates within the same month. The dates should also be consecutive, starting from the beginning of the month and going to the end. If the user enters start/end dates the cover more than one month, the data will all be recorded in a csv file named after the start date.
    If the user has to run the program multiple times for one month, the dates must follow each other consecutively. In other words, the program will not organize the csv files if the user has it run dates non-sequentially. Likewise, do not try to speed things up by splitting a month in hald and running the program twice at the same time for the same month. This will have give a csv non-sequential dates folded into each other, which will not behoove the user.
    You can, however, run the program on different months at the same time. So, you may, for example, have dates of 01/01/2015 - 01/31/2015 and 02/01/2015 - 02/28/2015 running at the same time. You cannot run 01/01/2015 - 01/15/2015 and 01/16/2015 at the same time, nor can you run 01/01/2015 - 02/28/2015


Other notes
    It takes about 48 hours for one whole month to have it's PFe counted. I advise runnig multiple months simultaneously.
    You must create a conda environment that can has all of the module below.
    If you encounter surprises, feel free to reach out at mcbarzal@terpmail.umd.edu. I may know whats going on.


Program outline
    This program reads in fits files from the hmi.Ic_noLimbDark_720s data series, accessed through heliocloud's aws S3 server. 
    It looks at the five hmi 12min (720s) averaged samples for each hour (i.e., for Sept 3, 2010, the first hour (midnight) has files from: [2010.09.03_00:00:00_TAI], [...00:12:00...], [...00:24:00...], [...00:36:00...], [...00:48:00...]). 
    The program checks to make sure all five files exist and do not have a quality flag (see jsoc.standford.edu). If all five files exist, they are derotated (because the sun rotates) to the 24min file and averaged together. These .fits are then saved as Avg_Map_date (program creates all directories for output files).
    The average map then has a mask applied to only the north pole and the south pole. These, too, are saved as north and south .fits files
    Photutils image segmentation is then applied to enhance and count the polar faculae in each mask. The scaling for the enchancement is borrowed from Munoz-Jaramillo et. al (2012) and Zhang (2010)
    The counts are appended into a csv file according the month. The csv also records any errors in data, so that every month has a row of data for every hour of every day, regardless of the data being useful.
    This results in months with 30 days having a csv of 720 rows and 31 days with 744 rows. If the month's csv does not have the right number of rows, the program will not be able to continue to the next part.
    
    A csv file is generated with the 720 and 744 data points smoothed to 10 north and south rows of data, with standard deviation.
    Plots are generated for each month, plotting Hourly Count vs Time, Smoothed Count vs Time, and Smoothed & Hourly vs Time.
    
'''



###############################################################################################
###############################################################################################
###############################################################################################


### Import modules
## Cell 1

import os
import io
import csv
import glob
import time
import warnings

import boto3
from botocore import UNSIGNED
from botocore.config import Config
import s3fs

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import datetime
from datetime import datetime as dt
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
import pandas as pd

import astropy.io.fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.fits import CompImageHDU

from photutils.segmentation import detect_sources,deblend_sources
from photutils.segmentation import SourceCatalog

import sunpy.map
import sunpy.data.sample
from sunpy.coordinates import frames, propagate_with_solar_surface, Helioprojective
from sunpy.map.maputils import all_coordinates_from_map, _verify_coordinate_helioprojective, solar_angular_radius


###############################################################################################
###############################################################################################
###############################################################################################


## Program intro  
print('\n\nSDO Data starts 05/20/2010. Do not enter an earlier date.\n')

## Gets user entered start and end dates, verifies they are in valid format and range
## Valid format is MM/DD/YYYY, range is from 05/20/2010 to present
while True: # Input start date, checks validity
    start_date = input('\nEnter start date, format MM/DD/YYYY\n')
    if len(start_date) == 10 and int(start_date[:2]) <= 12 and int(start_date[:2]) >= 1 and int(start_date[3:5]) >= 1 and int(start_date[3:5]) <= 31 and int(start_date[6:10]) >= 2010:
        break    # Good date!
    else:
        print('\nError: Invalid date.')     # Bad date! (we've all been there)
        
        
while True: # Input end date, checks validity
    end_date = input('\n\nEnter end date, format MM/DD/YYYY\n')
    if len(end_date) == 10 and int(end_date[:2]) <= 12 and int(end_date[:2]) >= 1 and int(end_date[3:5]) >= 1 and int(end_date[3:5]) <= 31 and int(end_date[6:10]) >= 2010:
        break     # Good date!
    else:
        print('\nError: invalid date.')     # Bad date!


## Double checks with user that entered dates are correct. 
## Cause if they arent, there will be a handful of directories
## to remove... not to mention, program takes a long time
while True:
    confirm_date = input('\n\nYou entered: start {0} end {1} \n Are these correct? y/n: ' .format(start_date,end_date))
    if confirm_date != 'n':
        break
    else: # If user says date is wrong, user may re-enter dates.
        
        while True: # Input start date, checks validity
            start_date = input('\nEnter start date, format MM/DD/YYYY\n')
            if len(start_date) == 10 and int(start_date[:2]) <= 12 and int(start_date[:2]) >= 1 and int(start_date[3:5]) >= 1 and int(start_date[3:5]) <= 31 and int(start_date[6:10]) >= 2010:
                break    # Good date!
            else:
                print('\nError: Invalid date.')     # Bad date! (we've all been there)

        while True: # Input end date, checks validity
            end_date = input('\n\nEnter end date, format MM/DD/YYYY\n')
            if len(end_date) == 10 and int(end_date[:2]) <= 12 and int(end_date[:2]) >= 1 and int(end_date[3:5]) >= 1 and int(end_date[3:5]) <= 31 and int(end_date[6:10]) >= 2010:
                break     # Good date!
            else:
                print('\nError: invalid date.')     # Bad date!
                
path = '/home/jovyan/efs/jrclark'
print('Ouput directory assumed /home/jovyan/efs/jrclark. See lines 145-149 if not.')

######
###### If you would like to output files to a different directory, 
###### you may edit 'path' at line 141 with your desired directory. 
###### If you would like have the program prompt you for an output
###### directory each time the program is run, comment out lines 141 
###### and 142, then remove the triple quotes on lines 152 and 160
######

'''
## Input a directory for program to poop files/figures into
while True:     # Check if user entered path exists
    path = input('\nEnter full path to directory in which you wish all generated files be placed.\n')
    if os.path.exists(path) == True: 
        break
    else:
        print('\nError: No such directory. Enter a path to an existing directory.')
'''
    
## If user entered their path with a "/" at the end, this
#  removes it so that I can enter my own slash
## Starting point of 'no slash' will make it easier to 
#  create individual bins/directories for different files to end up.
while path[-1] == '/':
    path = path[:-1]

    
## Ask if we're generating figures today. # variable used at _______________________
You_want_figures_with_that = input('\nYou want figures with that? y/n: ')
    
 
 ## Checks to see if user has used this program before and/or already has directories made up for all files
# If user does not yet have an organized directory, then I commit housekeeping.

if not os.path.exists(path + '/Saved FITS'):
    os.makedirs(path + '/Saved FITS') # Dir. for all FITS
if not os.path.exists(path + '/Saved FITS/Northern PFe Masks'):
    os.makedirs(path + '/Saved FITS/Northern PFe Masks') # Dir. for North Pole
if not os.path.exists(path + '/Saved FITS/Southern PFe Masks'):
    os.makedirs(path + '/Saved FITS/Southern PFe Masks') # Dir. for South Pole
if not os.path.exists(path + '/Saved FITS/Avg Maps'):
    os.makedirs(path + '/Saved FITS/Avg Maps') # Dir. for Avg

# Make directory for Year/Month fits files
yearmonth_dir = start_date[6:10] + '_' + start_date[:2]
if not os.path.exists(path + '/Saved FITS/Northern PFe Masks/' + yearmonth_dir):
    os.makedirs(path + '/Saved FITS/Northern PFe Masks/' + yearmonth_dir)
if not os.path.exists(path + '/Saved FITS/Southern PFe Masks/' + yearmonth_dir):
    os.makedirs(path + '/Saved FITS/Southern PFe Masks/' + yearmonth_dir) 
if not os.path.exists(path + '/Saved FITS/Avg Maps/' + yearmonth_dir):
    os.makedirs(path + '/Saved FITS/Avg Maps/' + yearmonth_dir) 

# Make Directory for saved .csv files
if not os.path.exists(path + '/PFe Data'):
    os.makedirs(path + '/PFe Data')
if not os.path.exists(path + '/PFe Data/Hourly'):
    os.makedirs(path + '/PFe Data/Hourly')
if not os.path.exists(path + '/PFe Data/Smoothed'):
    os.makedirs(path + '/PFe Data/Smoothed')  

    
###############################################################################################
###############################################################################################
###############################################################################################
    

## loading animation
'''
animation = [
"[        ]",
"[=       ]",
"[===     ]",
"[====    ]",
"[=====   ]",
"[======  ]",
"[======= ]",
"[========]",
"[ =======]",
"[  ======]",
"[   =====]",
"[    ====]",
"[     ===]",
"[      ==]",
"[       =]",
"[        ]",
"[        ]"
]
'''
animation = [
"[        ]",
"[☀       ]",
"[ ☀      ]",
"[  ☀     ]",
"[   ☀    ]",
"[    ☀   ]",
"[     ☀  ]",
"[      ☀ ]",
"[       ☀]",
"[      ☀ ]",
"[     ☀  ]",
"[    ☀   ]",
"[   ☀    ]",
"[  ☀     ]",
"[ ☀      ]",
"[☀       ]",
"[        ]",
"[☽       ]",
"[ ☽      ]",
"[  ☽     ]",
"[   ☽    ]",
"[    ☽   ]",
"[     ☽  ]",
"[      ☽ ]",
"[       ☽]",
"[      ☽ ]",
"[     ☽  ]",
"[    ☽   ]",
"[   ☽    ]",
"[  ☽     ]",
"[ ☽      ]",
"[☽       ]"
]

# animates
bling =0
def loading():
    global bling
    print(animation[bling % len(animation)], end='\r')
    time.sleep(.1)
    bling += 1
    return


###############################################################################################
###############################################################################################
###############################################################################################

# i don't know how else to call the loading animation, and it's not important enough to me to figure out right now
loading()
loading()
loading()
loading()
loading()
loading()
loading()
loading()

## Cell 2: Gather Data Keys and Group by Hour for Averaging
# Defines function to make list of grouped keys to access from data manifest files
# Data here is accessed through aws and S3
def get_grouped_keys(start,end):

    """ 
        Given a start and end date, will produce the list of keys for all Ic_noLimbDark_720s_hc files
        in that date range grouped by hour.

        start:Start date in the format 'M/D/YYYY' or 'MM/DD/YYYY'

        end: End date in in the format 'M/D/YYYY' or 'MM/DD/YYYY'

        return: A list of lists containing six objects. The first object is the time the FITS files
        in that list will eventually be rotated to. The other five items are the keys of all data within
        a specific hour.

        example: 

        keys_by_hour = get_grouped_keys('03/6/2022','3/10/2022'
        keys_by_hour

        [['2022-03-06T00:22:33.300',
        'sdo/hmi/20220306/Ic_noLimbDark/sdo_hmi_h2_20220306T000000_Ic_noLimbDark_v1.fits',
        'sdo/hmi/20220306/Ic_noLimbDark/sdo_hmi_h2_20220306T001200_Ic_noLimbDark_v1.fits',
        'sdo/hmi/20220306/Ic_noLimbDark/sdo_hmi_h2_20220306T002400_Ic_noLimbDark_v1.fits',
        'sdo/hmi/20220306/Ic_noLimbDark/sdo_hmi_h2_20220306T003600_Ic_noLimbDark_v1.fits',
        'sdo/hmi/20220306/Ic_noLimbDark/sdo_hmi_h2_20220306T004800_Ic_noLimbDark_v1.fits']...]

    """

    # Get list of dates in the correct format
    datelist = pd.date_range(start=start, end=end)
    dl = []
    for x in datelist:
        dl.append(str(x)[0:4]+str(x)[5:7]+str(x)[8:10])


    # Create list of manifest files for the given days
    manifest_ls = []
    for d in dl:
        manifest_ls.append('s3://gov-nasa-hdrl-data1/sdo/lists/hmi/' + str(d))

    # Make empty list
    keys_grouped = []
    # Individual files have hours indicated by T00 -> T23
    hours = ['T00', 'T01', 'T02', 'T03', 'T04', 'T05', 'T06', 'T07', 'T08', 'T09', 'T10', 'T11', 'T12', 'T13', 'T14', 'T15', 'T16', 'T17', 'T18', 'T19', 'T20', 'T21', 'T22', 'T23',]

    # Read manifest file and sort out only Ic_noLimbDark_720s_hc files
    for f in manifest_ls:
        day = str(f[39:43])+ '-' + str(f[43:45]) + '-' + str(f[45:47])
        Y = int(f[39:43])

        if Y >= 2013: 
            df_All = pd.read_csv(f)
            df_All.columns =['Data_Keys', 'Staging', 'Data Series', 'Date',]

            df_limbdark = df_All[df_All['Data Series'] == ' hmi.Ic_noLimbDark_720s_hc']

        if Y < 2013: # Files before 2013 were written different
            df_All = pd.read_csv(f)
            df_All.columns = ['Data_Keys', 'Staging', 'Date',]


            df_limbdark = df_All[df_All['Data_Keys'].str.contains('Ic_noLimbDark')] 

        for i in hours:
            keys_in_hour = []
            keys_in_hour.append(day+i+':22:33.300')
            for x in df_limbdark['Data_Keys']:
                if x[50:53] == i:
                    keys_in_hour.append(x)
            keys_grouped.append(keys_in_hour)

    return keys_grouped


loading()
###############################################################################################
###############################################################################################
###############################################################################################
loading()


## Cell 3

keys_by_hour = get_grouped_keys(start_date,end_date)


### Counting Algorithm

warnings.filterwarnings('ignore') # So that we don't see I/Iavg warning, SunPy Error

date_object = dt.strptime(start_date, "%m/%d/%Y") # Dates
csv_filename = date_object.strftime("/%b_%Y_Counts.csv") # Dates in good format, as filename

# If there is already a csv for this month, program will append it. Else: will create it anew
if not os.path.exists(path + '/PFe Data/Hourly' + csv_filename): 
    dfnum = pd.DataFrame(columns = ['Datetime', 'Day', 'N_Maps', 'North Count', 'South Count', 'Error'])
    dfnum.to_csv(path + '/PFe Data/Hourly' + csv_filename)


    
#Connect to S3, data stored here
mybucket='gov-nasa-hdrl-data1'
s3_res = boto3.resource('s3')
s3_bucket = s3_res.Bucket(mybucket)

s3c = boto3.client('s3', config=Config(signature_version=UNSIGNED)) # Version must be Unsigned, otherwise no cigar

#Define masks to be used later, limits data reading to >70deg latitudes (polar region)
def South_mask(Image_name):
    """ 
    Creates a mask covering pixels off the solar disk, the edge of the limb, and areas outside
    of the southern polar region (defined as anything above -70 degrees latitude)

    Image_name = SunPy Map

    Returns: an array the same size as the data showing ONLY the South Pole on the Sun

    """
    hpc_coords = all_coordinates_from_map(Image_name)
    _verify_coordinate_helioprojective(hpc_coords)
    all_hgs = hpc_coords.transform_to("heliographic_stonyhurst")

    mask = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) > (solar_angular_radius(hpc_coords)-3*u.arcsec)
    mask |= np.logical_or(all_hgs.lat >= -70 * u.deg, all_hgs.lat <= -90 * u.deg) 

    return mask


def North_mask(Image_name):
    """ 
    Creates a mask covering pixels off the solar disk, the edge of the limb, and areas outside
    of the northern polar region (defined as anything below 70 degrees latitude)

    Image_name = SunPy Map

    Returns: an array the same size as the data showing ONLY the South Pole on the Sun        
    """
    hpc_coords = all_coordinates_from_map(Image_name)
    _verify_coordinate_helioprojective(hpc_coords)
    all_hgs = hpc_coords.transform_to("heliographic_stonyhurst")

    mask = np.sqrt(hpc_coords.Tx ** 2 + hpc_coords.Ty ** 2) > (solar_angular_radius(hpc_coords)-3*u.arcsec)
    mask |= np.logical_or(all_hgs.lat >= 90 * u.deg, all_hgs.lat <= 70 * u.deg) 

    return mask


loading()
###############################################################################################
###############################################################################################
###############################################################################################
loading()


#This is where the code starts
for x in keys_by_hour:
    loading()
    time24 = x[0]

    #Check image 24 exists in hour
    check24 = False
    file24 = 'none'
    for i in x:
        loading()
        if i[53:55] == '24':
            check24 = True
            file24 = i

    #If Image 24 exists continue on to open fits files
    #If not give error: 'No Map24'
    if check24 == True:
        loading()
        MapList = [] #reprojected maps
        G_list = [] #All good quality maps

        try:
            loading()
            #Make Im24 which will be used to reproject other images later
            fobj24 = s3c.get_object(Bucket=mybucket,Key= file24)
            rawdata24 = fobj24['Body'].read()
            bdata24 = io.BytesIO(rawdata24)
            loading()
            hdul24 = astropy.io.fits.open(bdata24,memmap=False)
            if hdul24[1].header['QUALITY'] == 0:
                loading()
                Im24 = sunpy.map.Map(hdul24[1].data, hdul24[1].header)
                G_list.append(Im24.data)
                MapList.append(Im24.data)


            #Make other maps if quality is good
            for k in x:
                loading()
                if k[53:55] == '00':
                    loading()
                    fobj0 = s3c.get_object(Bucket=mybucket,Key= k)
                    rawdata0 = fobj0['Body'].read()
                    bdata0 = io.BytesIO(rawdata0)
                    loading()
                    hdul0 = astropy.io.fits.open(bdata0,memmap=False)
                    if hdul0[1].header['QUALITY'] == 0:
                        loading()
                        Im0 = sunpy.map.Map(hdul0[1].data, hdul0[1].header)
                        G_list.append(Im0)                        

                if k[53:55] == '12':
                    loading()
                    fobj12 = s3c.get_object(Bucket=mybucket,Key= k)
                    rawdata12 = fobj12['Body'].read()
                    bdata12 = io.BytesIO(rawdata12)
                    loading()
                    hdul12 = astropy.io.fits.open(bdata12,memmap=False)
                    if hdul12[1].header['QUALITY'] == 0:
                        loading()
                        Im12 = sunpy.map.Map(hdul12[1].data, hdul12[1].header)
                        G_list.append(Im12)

                if k[53:55] == '36':
                    loading()
                    fobj36 = s3c.get_object(Bucket=mybucket,Key= k)
                    rawdata36 = fobj36['Body'].read()
                    bdata36 = io.BytesIO(rawdata36)
                    loading()
                    hdul36 = astropy.io.fits.open(bdata36,memmap=False)
                    if hdul36[1].header['QUALITY'] == 0:
                        loading()
                        Im36 = sunpy.map.Map(hdul36[1].data, hdul36[1].header)
                        G_list.append(Im36)

                if k[53:55] == '48':
                    loading()
                    fobj48 = s3c.get_object(Bucket=mybucket,Key= k)
                    rawdata48 = fobj48['Body'].read()
                    bdata48 = io.BytesIO(rawdata48)
                    loading()
                    hdul48 = astropy.io.fits.open(bdata48,memmap=False)
                    if hdul48[1].header['QUALITY'] == 0:
                        loading()
                        Im48 = sunpy.map.Map(hdul48[1].data, hdul48[1].header)
                        G_list.append(Im48)

            #Check there are at least 5 good maps (derotate and average)
            #If not give error: 'N_Maps < 5'
            if len(G_list) >=5:
                loading()
                with propagate_with_solar_surface():
                    ProjIm0 = Im0.reproject_to(Im24.wcs)
                    MapList.append(ProjIm0.data)
                    ProjIm12 = Im12.reproject_to(Im24.wcs)
                    MapList.append(ProjIm12.data)
                    ProjIm36 = Im36.reproject_to(Im24.wcs)
                    MapList.append(ProjIm36.data)
                    ProjIm48 = Im48.reproject_to(Im24.wcs)
                    MapList.append(ProjIm48.data)           


                Output = 0
                for arr in MapList:
                    
                    Output = np.add(Output, arr)

                loading()
                AvgData = Output/len(MapList)

                # Save Averaged Map for later use
                # Get rid of NaN values. Replacing them by -1 will not interfere with the segmentation down below, 
                # unless your original data have some negative values
                AvgData[np.isnan(AvgData)] = -1
                # Convert to Float 32 bits
                AvgData = AvgData.astype(np.float32)
                AvgMap = sunpy.map.Map(AvgData,Im24.wcs)

                # Write the fits file with RICE compression
                loading()
                Avg_fname = path + '/Saved FITS/Avg Maps/' + yearmonth_dir + '/Avg_Map_' + time24 + '.fits'  
                AvgMap.save(Avg_fname, hdu_type=CompImageHDU, overwrite=True)
                loading()

                ImScaled = (AvgMap**15)*100

                SMask_fname = path + '/Saved FITS/Southern PFe Masks/' + yearmonth_dir + '/PFe_South_Mask_'+ time24 + '.fits' 
                NMask_fname = path + '/Saved FITS/Northern PFe Masks/' + yearmonth_dir + '/PFe_North_Mask_'+ time24 + '.fits'

                #Run image segmentation for the Northern Polar Region
                try:
                    loading()
                    MaskN = North_mask(ImScaled)
                    segment_mapN = detect_sources(ImScaled.data, 160, npixels=2, mask =MaskN)
                    N_Count = len(segment_mapN.labels)
                    chduN = fits.CompImageHDU(data=segment_mapN.data.astype(np.uint16), header=hdul24[1].header, compression_type='RICE_1')
                    chduN.writeto(NMask_fname, overwrite=True)
                    loading()

                except:
                    loading()
                    N_Count = np.nan

                #Run image segmentation for the Southern Polar Region
                try:
                    loading()
                    MaskS = South_mask(ImScaled)
                    segment_mapS = detect_sources(ImScaled.data, 160, npixels=2, mask =MaskS)
                    S_Count = len(segment_mapS.labels)
                    chduS = fits.CompImageHDU(data=segment_mapS.data.astype(np.uint16), header=hdul24[1].header, compression_type='RICE_1')
                    chduS.writeto(SMask_fname, overwrite=True)
                    loading()

                except:
                    loading()
                    S_Count = np.nan   

                new_row = pd.DataFrame({'Datetime': [time24], 'Day': [time24[:10]], 'N_Maps': [len(G_list)] , 'North Count': [N_Count], 'South Count': [S_Count], 'Error': ['None'] }) 
                new_row.to_csv(path + '/PFe Data/Hourly' + csv_filename, mode='a', index=True, header=False)
                loading()

            else:
                new_row = pd.DataFrame({'Datetime' : [time24], 'Day': [time24[:10]], 'N_Maps': [len(G_list)] , 'North Count': [np.nan], 'South Count': [np.nan], 'Error' : ['N_Maps < 5']})
                new_row.to_csv(path + '/PFe Data/Hourly' + csv_filename, mode='a', index=True, header=False)
                loading()
        except:
            new_row = pd.DataFrame({'Datetime' : [time24] , 'Day': [time24[:10]], 'N_Maps': [np.nan] , 'North Count': [np.nan], 'South Count': [np.nan], 'Error' : ['keys_DNE']})
            new_row.to_csv(path + '/PFe Data/Hourly' + csv_filename, mode='a', index=True, header=False)  
            loading()

    else:
        new_row = pd.DataFrame({'Datetime' : [time24] , 'Day': [time24[:10]], 'N_Maps': [np.nan] , 'North Count': [np.nan], 'South Count': [np.nan], 'Error' : ['No Map24']})
        new_row.to_csv(path + '/PFe Data/Hourly' + csv_filename, mode='a', index=True, header=False)
        loading()


warnings.filterwarnings('default') # Warnings back online


loading()
loading()
###############################################################################################
###############################################################################################
###############################################################################################


#Smoothing Data

#Read in Data here
df = pd.read_csv(path + '/PFe Data/Hourly' + csv_filename)

csv_smoothed_filename = date_object.strftime("/%b_%Y_Smoothed.csv")


dfstats = pd.DataFrame(columns = ['Date', 'N_North' , 'North Mean', 'North Stdv', 'N_South', 'South Mean', 'South Stdv'])


#Get all unique days in list
DayList = df['Day'].unique()

#Group days by 3 (we are averaging the data every 3 days)
daylist3 = []
for i in range(0, len(DayList), 3):
    chunk = DayList[i:i+3]
    daylist3.append(chunk)

#Average the Data
for x in daylist3:
    if len(x) == 3:
        df1 = df[df['Day'] == x[0]]
        df2 = df[df['Day'] == x[1]]
        df3 = df[df['Day'] == x[2]]

        dailydf = pd.concat([df1, df2, df3], ignore_index=True)

        Date = x[1]

        N_North = dailydf['North Count'].count()
        N_South = dailydf['South Count'].count()


        N_Mean = np.nanmean(dailydf['North Count'])
        N_stdv = np.nanstd(dailydf['North Count'])

        S_Mean = np.nanmean(dailydf['South Count'])
        S_stdv = np.nanstd(dailydf['South Count'])


        new_row = pd.Series({'Date': Date, 'N_North':N_North , 'North Mean': N_Mean,
                             'North Stdv': N_stdv, 'N_South': N_South, 'South Mean': S_Mean,
                             'South Stdv':S_stdv})

        dfstats = pd.concat([dfstats, new_row.to_frame().T], ignore_index=True) 

dfstats.to_csv(path + '/PFe Data/Smoothed' + csv_smoothed_filename)
loading()
loading()

##################################################################################################
##################################################################################################
##################################################################################################


#Do ya figure?

Year = date_object.strftime('%Y')
Year = int(Year)
Month = date_object.strftime('%B')

if You_want_figures_with_that != 'n':

    #Set up Date range for x-axis

    if Month == 'January':
        #Hourly
        date1 = datetime.datetime(Year, 1, 1, 0 ,22)
        date2 = datetime.datetime(Year, 2, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 1, 2, 0 ,22)
        date4 = datetime.datetime(Year, 1, 31, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Jan. '+ str(Year)
        save_date = 'Jan_'+ str(Year)

    if Month == 'February':
        if Year == 2012 or Year == 2016 or Year == 2020 or Year == 2024 or Year == 2028 or Year == 2032 or Year == 2036 or Year == 2040: #if yall using this past 2040, you wild.
            #Hourly
            date1 = datetime.datetime(Year, 2, 1, 0 ,22)
            date2 = datetime.datetime(Year, 3, 1, 0, 22)
            delta = datetime.timedelta(minutes=60)
            dates = drange(date1, date2, delta)

            #Smoothed
            date3 = datetime.datetime(Year, 2, 2, 0 ,22)
            date4 = datetime.datetime(Year, 2, 29, 0, 22)
            delta1 = datetime.timedelta(hours=72)
            dates1 = drange(date3, date4, delta1)

            title_date = 'Feb. '+ str(Year)
            save_date = 'Feb_'+ str(Year)
        else:
            #Hourly
            date1 = datetime.datetime(Year, 2, 1, 0 ,22)
            date2 = datetime.datetime(Year, 3, 1, 0, 22)
            delta = datetime.timedelta(minutes=60)
            dates = drange(date1, date2, delta)

            #Smoothed
            date3 = datetime.datetime(Year, 2, 2, 0 ,22)
            date4 = datetime.datetime(Year, 2, 28, 0, 22)
            delta1 = datetime.timedelta(hours=72)
            dates1 = drange(date3, date4, delta1)

            title_date = 'Feb. '+ str(Year)
            save_date = 'Feb_'+ str(Year)

    if Month == 'March':
        #Hourly
        date1 = datetime.datetime(Year, 3, 1, 0 ,22)
        date2 = datetime.datetime(Year, 4, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 3, 2, 0 ,22)
        date4 = datetime.datetime(Year, 3, 31, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Mar. '+ str(Year)
        save_date = 'Mar_'+ str(Year)

    if Month == 'April':
        #Hourly
        date1 = datetime.datetime(Year, 4, 1, 0, 22)
        date2 = datetime.datetime(Year, 5, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 4, 2, 0, 22)
        date4= datetime.datetime(Year, 4, 30, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Apr. '+ str(Year)
        save_date = 'Apr_'+ str(Year)

    if Month == 'May':
        #Hourly
        date1 = datetime.datetime(Year, 5, 1, 0 ,22)
        date2 = datetime.datetime(Year, 6, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 5, 2, 0 ,22)
        date4 = datetime.datetime(Year, 5, 31, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'May '+ str(Year)
        save_date = 'May_'+ str(Year)

    if Month == 'June':
        #Hourly
        date1 = datetime.datetime(Year, 6, 1, 0 ,22)
        date2 = datetime.datetime(Year, 7, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 6, 2, 0 ,22)
        date4 = datetime.datetime(Year, 6, 30, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Jun. '+ str(Year)
        save_date = 'Jun_'+ str(Year)

    if Month == 'July':
        #Hourly
        date1 = datetime.datetime(Year, 7, 1, 0 ,22)
        date2 = datetime.datetime(Year, 8, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 7, 2, 0 ,22)
        date4 = datetime.datetime(Year, 7, 31, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Jul. '+ str(Year)
        save_date = 'Jul_'+ str(Year)

    if Month == 'August':
        #Hourly
        date1 = datetime.datetime(Year, 8, 1, 0 ,22)
        date2 = datetime.datetime(Year, 9, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 8, 2, 0 ,22)
        date4 = datetime.datetime(Year, 8, 31, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Aug. '+ str(Year)
        save_date = 'Aug_'+ str(Year)

    if Month == 'September':
        #Hourly
        date1 = datetime.datetime(Year, 9, 1, 0 ,22)
        date2 = datetime.datetime(Year, 10, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 9, 2, 0 ,22)
        date4 = datetime.datetime(Year, 9, 30, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Sep. '+ str(Year)
        save_date = 'Sep_'+ str(Year)

    if Month == 'October':
        #Hourly
        date1 = datetime.datetime(Year, 10, 1, 0 ,22)
        date2 = datetime.datetime(Year, 11, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 10, 2, 0 ,22)
        date4 = datetime.datetime(Year, 10, 31, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Oct. '+ str(Year)
        save_date = 'Oct_'+ str(Year)

    if Month == 'November':
        #Hourly
        date1 = datetime.datetime(Year, 11, 1, 0 ,22)
        date2 = datetime.datetime(Year, 12, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 11, 2, 0 ,22)
        date4 = datetime.datetime(Year, 11, 30, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Nov. '+ str(Year)
        save_date = 'Nov_'+ str(Year)

    if Month == 'December':
        #Hourly
        date1 = datetime.datetime(Year, 12, 1, 0 ,22)
        date2 = datetime.datetime(Year + 1, 1, 1, 0, 22)
        delta = datetime.timedelta(minutes=60)
        dates = drange(date1, date2, delta)

        #Smoothed
        date3 = datetime.datetime(Year, 12, 2, 0 ,22)
        date4 = datetime.datetime(Year, 12, 31, 0, 22)
        delta1 = datetime.timedelta(hours=72)
        dates1 = drange(date3, date4, delta1)

        title_date = 'Dec. '+ str(Year)
        save_date = 'Dec_'+ str(Year)


    loading()


    if not os.path.exists(path + '/Figures'):
        os.makedirs(path + '/Figures')
    if not os.path.exists(path + '/Figures/Hourly'):
        os.makedirs(path + '/Figures/Hourly')
    if not os.path.exists(path + '/Figures/Smoothed'):
        os.makedirs(path + '/Figures/Smoothed')
    if not os.path.exists(path + '/Figures/Smoothed and Hourly'):
        os.makedirs(path + '/Figures/Smoothed and Hourly')

        
    #Plot of PFe Counts over time

    fig, ax = plt.subplots()
    ax.plot(dates, df['North Count'], color='#19ceef',linestyle = '-', marker='D', markersize=3,label = 'PFe Count North')
    ax.plot(dates, df['South Count'], color='#2f2da4',linestyle = '-', marker='s', markersize=3,label = 'PFe Count South')

    ax.xaxis.set_major_locator(DayLocator(range(1,30,5)))
    ax.xaxis.set_minor_locator(DayLocator(range(1,30)))
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    fig.autofmt_xdate()

    ax.legend(bbox_to_anchor = (1, 1))

    ax.set_ylabel('PFe Count', color='k')
    ax.set_xlabel('Time (UTC)')
    plt.title('PFe Count in North and South [' + title_date + ']')

    plt.savefig(path + '/Figures/Hourly/' + save_date + '_Data.png',bbox_inches = 'tight')

    loading()

    
    #Smoothed Data
 
    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot(dates1, dfstats['North Mean'], color='#19ceef',linestyle = '-', marker='D', markersize=4,
            label = 'Northern PFe Count')
    ax.plot(dates1, dfstats['South Mean'], color='#2f2da4',linestyle = '-', marker='s', markersize=4,
            label = 'South PFe Count')

    ax.fill_between(dates1, np.array(dfstats['North Mean'] - dfstats['North Stdv'], dtype=float), np.array(dfstats['North Mean'] + dfstats['North Stdv'], dtype=float),
                   color='#19ceef', alpha = 0.2, label = 'North stdv')

    ax.fill_between(dates1, np.array(dfstats['South Mean'] - dfstats['South Stdv'], dtype=float), np.array(dfstats['South Mean'] + dfstats['South Stdv'], dtype=float),
                   color='#2f2da4', alpha = 0.2, label = 'South stdv')

    ax.xaxis.set_major_locator(DayLocator(range(1,30,5)))
    ax.xaxis.set_minor_locator(DayLocator(range(1,30)))
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    fig.autofmt_xdate()

    ax.legend(bbox_to_anchor = (1,1),loc='upper left')

    ax.set_ylabel('PFe Count', color='k')
    ax.set_xlabel('Time')
    plt.title('Smoothed PFe Counts [' + title_date + ']')

    plt.savefig(path + '/Figures/Smoothed/' + save_date + '_Smoothed.png',bbox_inches = 'tight')

    loading()


    #Original and Smoothed

    fig, ax = plt.subplots(figsize=(7,5))
    ax.plot(dates, df['North Count'], color='#19ceef',linestyle = ' ', marker='D', markersize=4,label = 'PFe Count North')
    ax.plot(dates, df['South Count'], color='#2f2da4',linestyle = ' ', marker='s', markersize=4,label = 'PFe Count South')

    ax.fill_between(dates1, np.array(dfstats['North Mean'] - dfstats['North Stdv'], dtype=float), np.array(dfstats['North Mean'] + dfstats['North Stdv'], dtype=float),
                   color='#19ceef', alpha = 0.2, label = 'North stdv')

    ax.fill_between(dates1, np.array(dfstats['South Mean'] - dfstats['South Stdv'], dtype=float), np.array(dfstats['South Mean'] + dfstats['South Stdv'], dtype=float),
                   color='#2f2da4', alpha = 0.3, label = 'South stdv')

    ax.xaxis.set_major_locator(DayLocator(range(1,30,5)))
    ax.xaxis.set_minor_locator(DayLocator(range(1,30)))
    ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

    fig.autofmt_xdate()

    ax.legend(bbox_to_anchor = (1,1),loc='upper left')

    ax.set_ylabel('PFe Count', color='k')
    ax.set_xlabel('Time')
    plt.title('Smoothed and Hourly Data [' + title_date + ']')

    plt.savefig(path + '/Figures/Smoothed and Hourly/' + save_date + '_Smoothed_and_Hourly.png',bbox_inches = 'tight')

    loading()
    loading()
    loading()
    loading()

    print('\nDone')

else:
    loading()
    loading()
    loading()
    loading()

    print('\nDone')
