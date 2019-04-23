#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Read in raw TGI data
# 
# Version history:
#
# 02-Apr-2019  M. Peel       Started
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import pandas as pd
import scipy.fftpack
import matplotlib as mpl
from scipy import signal, optimize
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
from astropy.time import Time
from astropy.coordinates import Angle
from scipy.optimize import curve_fit
import tgi

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

# data = get_sci('testdata/CRAB-19-03-01-18-07-38-0000.sci2')
# plot_sci(data)
# print(data)
# exit()
# dt = np.dtype([('time', [('min', int), ('sec', int)]),('temp', float)])

# run = tgi.tgi(indir='/Volumes/proyectos/quijote2/tod/',outdir='/Users/mpeel/Documents/git/quijote')
run = tgi.tgi(indir='/Users/mpeel/Documents/git/quijote/testdata/',outdir='/Users/mpeel/Documents/git/quijote',pixelfileloc='/Users/mpeel/Documents/git/quijote/etc/qt2_masterfiles_new/ qt2_pixel_masterfile.',pixelposfileloc='/Users/mpeel/Documents/git/quijote/etc/tgi_fgi_horn_positions_table.txt')

# pixelinfo = run.get_pixel_info(2457679.0,16)
# print(pixelinfo)
# print(pixelinfo['fp'])
# exit()

# run.stack_maps_tod('CRAB',['CRAB-190410-2025'],pixelrange=[3,4,5,21,22,23,24],detrange=[0],phaserange=[0],plotlimit=0.001,numfiles=30,dopol=False)
# run.stack_maps_tod('HAZE',['HAZE-190411-0238','HAZE-190411-0626','HAZE-190412-0234','HAZE-190412-0622'],pixelrange=[3,4,5,21,22,23,24],detrange=[0],phaserange=[0],plotlimit=0.001,numfiles=50,dopol=False)
# exit()
# run.analyse_tod('CRAB-190311-1728',numfiles=1,pixelrange=[3,4,5,21,22,23,24],detrange=[0],dopol=True)
# exit()
run.analyse_tod('CRAB-190409-2029',dopol=True)
exit()
# tgi.analyse_tod('testdata/CRAB-190410-1529',numfiles=14)
# run.analyse_tod('CRAB-190410-2025',pixelrange=[23,24],detrange=[0],plotlimit=0.001,numfiles=50,dopol=True)
# run.analyse_tod('CRAB-190410-2025',pixelrange=[3,4,5,21,22,23,24],plotlimit=0.001,numfiles=1)
# for i in range(1,5):
# 	run.combine_sky_maps(['/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/skymap_4_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/skymap_5_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/skymap_6_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/skymap_24_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/skymap_25_1_'+str(i)+'.fits'],['/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/hitmap_4_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/hitmap_5_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/hitmap_6_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/hitmap_24_1_'+str(i)+'.fits','/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/hitmap_25_1_'+str(i)+'.fits'],'/Users/mpeel/Documents/git/quijote/CRAB-190410-2025_pol/combined_'+str(i),centralpos=(184,-5))

# run.analyse_tod('MOON-180919-1806',plotlimit=0.002,dopol=True)
run.analyse_tod('MOON-180921-2035',plotlimit=0.002,dopol=True)

if 1:
	name = 'MOON-180919-1806'
	centralpos=(299,-20)
	# name = 'MOON-180921-2035'
	# centralpos=(37,-44)
	for i in range(1,5):
		if i == 2 or i == 3:
			plotlimit = 0.0001
		else:
			plotlimit = 0.0
		for k in range(1,2):
			filelist = []
			hitlist = []
			for pix in range(1,31):
				filelist.append('/Users/mpeel/Documents/git/quijote/'+name+'/skymap_'+str(pix)+'_'+str(k)+'_'+str(i)+'.fits')
				hitlist.append('/Users/mpeel/Documents/git/quijote/'+name+'/hitmap_'+str(pix)+'_'+str(k)+'_'+str(i)+'.fits')
			print(filelist)
			run.combine_sky_maps(filelist,hitlist,'/Users/mpeel/Documents/git/quijote/'+name+'/combined2_'+str(i)+'_'+str(k),centralpos=centralpos,plotlimit=plotlimit)

#run.analyse_tod('MOON-180921-2035',plotlimit=0.001,numfiles=50,dopol=True,plottods=False)

exit()

# 2019-04-11
pixelrange=[3,4,5,21,22,23,24]

datasets = ['HAZE-190411-0238','CYGNUS-190411-0452','HAZE-190411-0626','CASS-190411-0840','CYGNUS-190411-0955','PERSEUS-190411-1134','PERSEUS-190411-1257','CASS-190411-1423','MOON-190411-1525','M42-190411-1631','PERSEUS-190411-1734','PERSEUS-190411-1858','CRAB-190411-2021']
#'HAZE-190412-0234','CYGNUS-190412-0448','HAZE-190412-0622',
datasets = ['CASS-190412-0836','CYGNUS-190412-0951','PERSEUS-190412-1131','PERSEUS-190412-1254','CASS-190412-1420']

#run.analyse_tod('DIP000-190411-2120',pixelrange=[3,4,5,21,22,23,24],plotlimit=0.001)
#run.analyse_tod('DIP000-190411-2143',pixelrange=[3,4,5,21,22,23,24],plotlimit=0.001)
for dataset in datasets:
	run.analyse_tod(dataset,pixelrange=pixelrange,plotlimit=0.001)
