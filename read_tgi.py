#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Run the TFGI python pipeline
# 
# Version history:
#
# 02-Apr-2019  M. Peel       Started
# 17-May-2019  M. Peel 		 Tidied up

import healpy as hp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tfgi

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

# Set this to where you have a copy of the data
# basedir = '/Volumes/proyectos/quijote2/'
basedir = '/Users/mpeel/Documents/git/quijote/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
outdir = '/Users/mpeel/Documents/git/quijote/output/'
# Start the class
run = tfgi.tfgi(outdir=outdir,\
	datadir=basedir+'tod/',\
	pixelfileloc=basedir+'etc/qt2_pixel_masterfile.',\
	pixelposfileloc=basedir+'etc/tgi_fgi_horn_positions_table.txt',\
	polcalfileloc=basedir+'etc/qt2_masterfiles_new/qt2_polcal.')

# frequencies = [150e9,220e9,11e9,30e9,40e9]
# diameters = [0.22,0.22,2.5,2.5,2.5]
# for i in range(0,len(frequencies)):
# 	print(run.calc_farfield(diameters[i],frequency=frequencies[i]))
# exit()



# run.analyse_tod('CYGNUS-190411-0452',numfiles=1)

# run.analyse_tod('MOON-190411-1525')
# exit()
# run.analyse_skydip('DIP000-190411-2120',dopol=True)
# exit()

# datasets = run.find_observations('MOON-19')
datasets = run.find_observations('DIP')
# datasets2 = run.find_observations('CRAB-1903')
# datasets = list(set(datasets1) | set(datasets2))
print(datasets)
for dataset in datasets:
	# You can set options for the reduction in the next line. The options and their defaults are:
	# pixelrange=range(0,31)  - set to an array, defaults to all pixels and the masterfile
	# detrange=range(0,4)     - set to an array, defaults to all detectors
	# phaserange=range(0,4)   - set to an array, defaults to all phase outputs
	# plotlimit=0.0           - lets you define a maximum value in some output plots
	# quiet=False             - set to true if you want the code to run quietly
	# dofft=False             - set to false to generate ffts and fit for f_knee
	# plottods=True           - set to false to not plot tods
	# plotmap=True            - set to false to not plot maps
	# dopol=False             - set to true to change from detector to polarised outputs
	# plotcombination=True    - set to false to disable creating a combined map
	# numfiles=50             - set to a lower number to only read in the first files of each observation
	run.analyse_skydip(dataset,dopol=False)
	# run.analyse_tod(dataset,plotlimit=0.001,dopol=False,plottods=False)

exit()

# Testing source features
# run.examine_source(['/Users/mpeel/Documents/git/quijote/MOON-190411-1525/skymap_24_1_1.fits'],['/Users/mpeel/Documents/git/quijote/MOON-190411-1525/hitmap_4_1_1.fits'],'test')

# Some test skydips
# run.analyse_skydip('DIP000-190411-2120',dopol=True,detrange=[0],numelbins=50)
# run.analyse_skydip('DIP000-190411-2143',dopol=True)

# exit()

# pixels = [[1]]
# channels = [[17]]
# files = [['PolCal_june_2018/June_22_2018/Pix1-18-06-22-10-25-28-0000.sci2','PolCal_june_2018/June_22_2018/Pix1-18-06-22-10-25-28-0001.sci2','PolCal_june_2018/June_22_2018/Pix1-18-06-22-10-25-28-0002.sci2','PolCal_june_2018/June_22_2018/Pix1-18-06-22-10-25-28-0003.sci2','PolCal_june_2018/June_22_2018/Pix1-18-06-22-10-25-28-0004.sci2','PolCal_june_2018/June_22_2018/Pix1-18-06-22-10-25-28-0005.sci2']]
# files = [['PolCal_june_2018/Pix1_cal-18-06-04-13-33-26-0000.sci2','PolCal_june_2018/Pix1_cal-18-06-04-13-33-26-0001.sci2','PolCal_june_2018/Pix1_cal-18-06-04-13-33-26-0002.sci2','PolCal_june_2018/Pix1_cal-18-06-04-13-33-26-0003.sci2','PolCal_june_2018/Pix1_cal-18-06-04-13-33-26-0004.sci2','PolCal_june_2018/Pix1_cal-18-06-04-13-33-26-0005.sci2']]

## April 2019 calibration data

# Run through the scientific data
pixels = [[5],[17],[23],[26],[41],[42],[63]]
fix_neg = [False,False,False,False,True,True,False]
channels = [[25],[23],[24],[22],[4],[6],[5]]
files = [['testdata/'+'pix5_cal-19-04-09-13-14-17-0000.sci2','testdata/'+'pix5_cal-19-04-09-13-14-17-0001.sci2'],['testdata/'+'pix17_cal-19-04-09-12-56-16-0000.sci2','testdata/'+'pix17_cal-19-04-09-12-56-16-0001.sci2'],['testdata/'+'pix23_cal-19-04-09-13-33-21-0000.sci2','testdata/'+'pix23_cal-19-04-09-13-33-21-0001.sci2'],['testdata/'+'pix26_cal-19-04-09-12-20-15-0000.sci2','testdata/'+'pix26_cal-19-04-09-12-20-15-0001.sci2'],['testdata/'+'Pix41_cal-19-04-10-11-31-14-0000.sci2','testdata/'+'Pix41_cal-19-04-10-11-31-14-0001.sci2'],['testdata/'+'Pix42_cal-19-04-10-12-07-14-0000.sci2','testdata/'+'Pix42_cal-19-04-10-12-07-14-0001.sci2'],['testdata/'+'Pix63_cal-19-04-10-11-48-14-0000.sci2','testdata/'+'Pix63_cal-19-04-10-11-48-14-0001.sci2']]
for i in range(0,len(pixels)):
	print('Pixel '+str(pixels[i]))
	data = run.get_sci(files[i],quiet=True)
	run.plot_sci(data,channels=channels[i],pixels=pixels[i],fix_neg=fix_neg[i],offset=True)
exit()

data = run.get_sci(['testdata/'+'pix5_cal-19-04-09-13-14-17-0000.sci2','testdata/'+'pix5_cal-19-04-09-13-14-17-0001.sci2'],quiet=True)
run.plot_sci(data,channels=[25],pixels=[5],fix_neg=False,offset=True)


# Run through the engineering data
pixels = [5,17,23,26,41,42,63]
files = ['pix5_cal-19-04-09-13-10-07-0000.eng2','pix17_cal-19-04-09-12-51-20-0000.eng2','pix23_cal-19-04-09-13-28-32-0000.eng2','pix26_cal-19-04-09-11-45-34-0000.eng2','Pix41_cal-19-04-10-11-26-25-0000.eng2','Pix42_cal-19-04-10-12-03-14-0000.eng2','Pix63_cal-19-04-10-11-43-50-0000.eng2']
for i in range(0,len(pixels)):
	print('Pixel '+str(pixels[i]))
	data = run.get_eng('testdata/'+files[i],quiet=True)
	run.plot_eng(data,pixel=pixels[i])
exit()
