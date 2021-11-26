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

## November 2021 calibration data

# Run through the scientific data
# pixels = [[5],[19],[23],[26],[41],[42],[63]]
# # channels = [[25],[23],[24],[21],[4],[6],[5]]
# files = [['testdata/pix41_polcal-21-11-22-12-25-54-0000.sci2','testdata/pix41_polcal-21-11-22-12-25-54-0001.sci2']]
# files = [['testdata/'+'pix5_cal-19-04-09-13-14-17-0000.sci2','testdata/'+'pix5_cal-19-04-09-13-14-17-0001.sci2'],['testdata/'+'pix17_cal-19-04-09-12-56-16-0000.sci2','testdata/'+'pix17_cal-19-04-09-12-56-16-0001.sci2'],['testdata/'+'pix23_cal-19-04-09-13-33-21-0000.sci2','testdata/'+'pix23_cal-19-04-09-13-33-21-0001.sci2'],['testdata/'+'pix26_cal-19-04-09-12-20-15-0000.sci2','testdata/'+'pix26_cal-19-04-09-12-20-15-0001.sci2'],['testdata/'+'Pix41_cal-19-04-10-11-31-14-0000.sci2','testdata/'+'Pix41_cal-19-04-10-11-31-14-0001.sci2'],['testdata/'+'Pix42_cal-19-04-10-12-07-14-0000.sci2','testdata/'+'Pix42_cal-19-04-10-12-07-14-0001.sci2'],['testdata/'+'Pix63_cal-19-04-10-11-48-14-0000.sci2','testdata/'+'Pix63_cal-19-04-10-11-48-14-0001.sci2']]

pixels = [[41]]
channels = [[4]]
# files = [['testdata/pix41_polcal-21-11-22-12-25-54-0000.sci2','testdata/pix41_polcal-21-11-22-12-25-54-0001.sci2']]
files = [['testdata/PIX_41_polcal-21-11-23-12-23-03-0000.sci2','testdata/PIX_41_polcal-21-11-23-12-23-03-0001.sci2']]

# pixels = [[63]]
# channels = [[5]]
# files = [['testdata/PIX_63_polcal-21-11-23-11-47-05-0000.sci2','testdata/PIX_63_polcal-21-11-23-11-47-05-0001.sci2']]
# pixels = [[42]]
# channels = [[6]]
# files = [['testdata/PIX_42_polcal-21-11-23-12-08-09-0000.sci2','testdata/PIX_42_polcal-21-11-23-12-08-09-0001.sci2']]
fix_neg = [False,False,False,False,True,True,False]
for i in range(0,len(pixels)):
	print('Pixel '+str(pixels[i]))
	data = run.get_sci(files[i],quiet=True)
	run.plot_sci(data,channels=channels[i],pixels=pixels[i],fix_neg=fix_neg[i],offset=True)

# exit()
#
# data = run.get_sci(['testdata/'+'pix5_cal-19-04-09-13-14-17-0000.sci2','testdata/'+'pix5_cal-19-04-09-13-14-17-0001.sci2'],quiet=True)
# run.plot_sci(data,channels=[25],pixels=[5],fix_neg=False,offset=True)


# Run through the engineering data
# pixels = [5,17,23,26,41,42,63]
# pixels = [41]
# files = ['pix41_polcal-21-11-22-12-21-19-0000.eng2','pix41_polcal-21-11-22-12-21-19-0001.eng2']
# files = ['PIX_41_polcal-21-11-23-12-20-01-0000.eng2','PIX_41_polcal-21-11-23-12-20-01-0001.eng2']
# pixels = [63]
# files = ['PIX_63_polcal-21-11-23-11-44-06-0000.eng2','PIX_63_polcal-21-11-23-11-44-06-0001.eng2']
# pixels = [42]
# files = ['PIX_42_polcal-21-11-23-12-04-01-0000.eng2','PIX_42_polcal-21-11-23-12-04-01-0001.eng2']
# pixels = [23]
# files = ['PIX_23_polcal-21-11-23-13-32-01-0000.eng2','PIX_23_polcal-21-11-23-13-32-01-0001.eng2']
# for i in range(0,len(pixels)):
# 	print('Pixel '+str(pixels[i]))
# 	data = run.get_eng('testdata/'+files[i],quiet=True)
# 	run.plot_eng(data,pixel=pixels[i])
# exit()
