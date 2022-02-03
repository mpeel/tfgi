#!/usr/bin/env python
# -*- coding: utf-8  -*-
import healpy as hp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tfgi

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

# Set this to where you have a copy of the data
# basedir = '/Volumes/proyectos/quijote2/'
basedir = '/Volumes/WD12TB/quijote2/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
outdir = '/Volumes/WD12TB/quijote2/output/'
# Start the class
run = tfgi.tfgi(outdir=outdir,\
	datadir=basedir+'tod/',\
	pixelfileloc=basedir+'etc/qt2_pixel_masterfile.',\
	pixelposfileloc=basedir+'etc/tgi_fgi_horn_positions_table.txt',\
	polcalfileloc=basedir+'etc/qt2_masterfiles_new/qt2_polcal.')

# January 2022 calibration data

# Run through the scientific data
pixels = [[5], [5], [19], [19], [23], [23], [26], [26]]
channels = [[25], [25], [23], [23], [24], [24], [21], [21]]
files = [['pix_5_polcal1-22-01-12-12-17-03-0000.sci2','pix_5_polcal1-22-01-12-12-17-03-0001.sci2'], ['pix_5_polcal2-22-01-12-13-40-02-0000.sci2', 'pix_5_polcal2-22-01-12-13-40-02-0001.sci2'], ['pix_19_polcal1-22-01-12-11-32-02-0000.sci2', 'pix_19_polcal1-22-01-12-11-32-02-0001.sci2'], ['pix_19_polcal2-22-01-12-13-26-00-0000.sci2', 'pix_19_polcal2-22-01-12-13-26-00-0001.sci2'], ['pix_23_polcal1-22-01-12-11-47-09-0000.sci2', 'pix_23_polcal1-22-01-12-11-47-09-0001.sci2'], ['pix_23_polcal2-22-01-12-12-57-03-0000.sci2', 'pix_23_polcal2-22-01-12-12-57-03-0001.sci2'], ['pix_26_polcal1-22-01-12-12-03-03-0000.sci2', 'pix_26_polcal1-22-01-12-12-03-03-0001.sci2'], ['pix_26_polcal2-22-01-12-13-12-02-0000.sci2', 'pix_26_polcal2-22-01-12-13-12-02-0001.sci2']]
prefix = ['polcal1','polcal2', 'polcal1','polcal2', 'polcal1','polcal2', 'polcal1','polcal2']
fix_neg = [False,False,False,False,False,False,False,False]
for i in range(0,len(pixels)):
	if i < 6:
		continue
	print('Pixel '+str(pixels[i]))
	data = run.get_sci(files[i],quiet=True,indir=basedir+'data/2022-01/')
	try:
		run.plot_sci(data,channels=channels[i],pixels=pixels[i],fix_neg=fix_neg[i],offset=True,prefix=prefix[i])
	except:
		print('There was an error here')
		pass


# Run through the engineering data
# pixels = [5,17,23,26,41,42,63]
# pixels = [41]
# files = ['pix41_polcal-21-11-22-12-21-19-0000.eng2','pix41_polcal-21-11-22-12-21-19-0001.eng2']
# files = ['PIX_41_polcal-21-11-23-12-20-01-0000.eng2','PIX_41_polcal-21-11-23-12-20-01-0001.eng2']
# pixels = [63]
# files = ['PIX_63_polcal-21-11-23-11-44-06-0000.eng2','PIX_63_polcal-21-11-23-11-44-06-0001.eng2']
# pixels = [42]
# files = ['PIX_42_polcal-21-11-23-12-04-01-0000.eng2','PIX_42_polcal-21-11-23-12-04-01-0001.eng2']
# pixels = [5]
# files = ['pix_5_polcal1-22-01-12-12-14-04-0000.eng2','pix_5_polcal1-22-01-12-12-14-04-0001.eng2']#,'pix_5_polcal1-22-01-12-12-14-04-0002.eng2']
# for i in range(0,len(pixels)):
# 	print('Pixel '+str(pixels[i]))
# 	data = run.get_eng(basedir+'data/2022-01/'+files[i],quiet=True)
# 	run.plot_eng(data,pixel=pixels[i])
# exit()
