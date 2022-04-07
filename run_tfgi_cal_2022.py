#!/usr/bin/env python
# -*- coding: utf-8  -*-
import healpy as hp
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tfgi

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000

## Set this to where you have a copy of the data (with 'data' and 'etc' subdirectories)
# basedir = '/Volumes/WD12TB/quijote2/'
# basedir = '/Users/mpeel/Documents/git/tfgi/'
basedir = '/Volumes/proyectos/quijote2/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
# outdir = '/Volumes/WD12TB/quijote2/output/'
outdir = '/Users/mpeel/Documents/git/tfgi/'

# Start the class
run = tfgi.tfgi(outdir=outdir,\
	datadir=basedir+'tod/',\
	pixelfileloc=basedir+'etc/qt2_pixel_masterfile.',\
	pixelposfileloc=basedir+'etc/tgi_fgi_horn_positions_table.txt',\
	polcalfileloc=basedir+'etc/qt2_masterfiles_new/qt2_polcal.')

# Run through the scientific data
# NOTE: if these don't run, make sure you have a ,\ at the end of each line - EXCEPT the last one, where it should end with a }
measurements = [\

## 2022-04-07
{'pixel':63,'channel':5,'files':['pixel63_polcal1-22-04-07-14-12-03-0000.sci2','pixel63_polcal1-22-04-07-14-12-03-0001.sci2'],'prefix':'2022-04-07_polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':41,'channel':4,'files':['pixel41_polcal1-22-04-07-13-57-58-0000.sci2','pixel41_polcal1-22-04-07-13-57-58-0001.sci2'],'prefix':'2022-04-07_polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':5,'channel':25,'files':['pixel5_polcal2-22-04-07-12-52-03-0000.sci2','pixel5_polcal2-22-04-07-12-52-03-0001.sci2'],'prefix':'2022-04-07_polcal2_new','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':5,'channel':25,'files':['pixel5_polcal2-22-04-07-11-59-59-0000.sci2','pixel5_polcal2-22-04-07-11-59-59-0001.sci2'],'prefix':'2022-04-07_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # NOTE: NEXT FILES MISNAMED, SAYS 26 ACTUALLY 19
# # {'pixel':19,'channel':23,'files':['pixel26_polcal2-22-04-07-11-40-59-0000.sci2','pixel26_polcal2-22-04-07-11-40-59-0001.sci2'],'prefix':'2022-04-07_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\ # BAD MEASUREMENT?
# {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-07-11-24-59-0000.sci2','pixel26_polcal2-22-04-07-11-24-59-0001.sci2'],'prefix':'2022-04-07_polcal2_new','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-07-10-42-59-0000.sci2','pixel26_polcal2-22-04-07-10-42-59-0001.sci2'],'prefix':'2022-04-07_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\ # BAD MEASUREMENT?
# {'pixel':23,'channel':24,'files':['pixel23_polcal2-22-04-07-09-48-58-0000.sci2','pixel23_polcal2-22-04-07-09-48-58-0001.sci2'],'prefix':'2022-04-07_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\

## 2022-04-06
# {'pixel':5,'channel':25,'files':['pixel5_polcal2-22-04-06-16-03-58-0000.sci2','pixel5_polcal2-22-04-06-16-03-58-0001.sci2'],'prefix':'2022-04-06_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # {'pixel':19,'channel':23,'files':['pixel19_polcal2-22-04-06-15-46-00-0000.sci2','pixel19_polcal2-22-04-06-15-46-00-0001.sci2'],'prefix':'2022-04-06_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\ # MEASUREMENT BAD???
# # {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-06-15-31-58-0000.sci2','pixel26_polcal2-22-04-06-15-31-58-0001.sci2'],'prefix':'2022-04-06_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\ # MEASUREMENT BAD??
# # {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-06-15-24-02-0000.sci2','pixel26_polcal2-22-04-06-15-24-02-0001.sci2'],'prefix':'2022-04-06_polcal2_bad','fix_neg':False,'indir':basedir+'data/2022-04/'},\ # MEASUREMENT BAD??
# {'pixel':23,'channel':24,'files':['pixel23_polcal2-22-04-06-15-06-58-0000.sci2','pixel23_polcal2-22-04-06-15-06-58-0001.sci2'],'prefix':'2022-04-06_polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':23,'channel':24,'files':['pixel23_polcal2-22-04-06-14-58-58-0000.sci2','pixel23_polcal2-22-04-06-14-58-58-0001.sci2'],'prefix':'2022-04-06_polcal2_bad','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':5,'channel':25,'files':['pixel5_polcal1-22-04-06-13-59-57-0000.sci2','pixel5_polcal1-22-04-06-13-59-57-0001.sci2'],'prefix':'2022-04-06_polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':19,'channel':23,'files':['pixel19_polcal1-22-04-06-13-40-56-0000.sci2','pixel19_polcal1-22-04-06-13-40-56-0001.sci2'],'prefix':'2022-04-06_polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':26,'channel':21,'files':['pixel26_polcal1-22-04-06-13-26-00-0000.sci2','pixel26_polcal1-22-04-06-13-26-00-0001.sci2'],'prefix':'2022-04-06_polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':23,'channel':24,'files':['pixel23_polcal1-22-04-06-12-08-00-0000.sci2','pixel23_polcal1-22-04-06-12-08-00-0001.sci2'],'prefix':'2022-04-06_polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'}

## 2022-03-11
# {'pixel':5,'channel':25,'files':['pixel_5_polcal1-22-03-11-13-55-19-0000.sci2','pixel_5_polcal1-22-03-11-13-55-19-0001.sci2'],'prefix':'2022-03-11_polcal1','fix_neg':False,'indir':basedir+'data/2022-03/'},\
# {'pixel':19,'channel':23,'files':['pixel_19_polcal1-22-03-11-13-41-11-0000.sci2', 'pixel_19_polcal1-22-03-11-13-41-11-0001.sci2'],'prefix':'2022-03-11_polcal1','fix_neg':False,'indir':basedir+'data/2022-03/'},\
# {'pixel':23,'channel':24,'files':['pixel_23_polcal1-22-03-11-13-11-08-0000.sci2', 'pixel_23_polcal1-22-03-11-13-11-08-0001.sci2'],'prefix':'2022-03-11_polcal1','fix_neg':False,'indir':basedir+'data/2022-03/'},\
## {'pixel':26,'channel':21,'files':['pixel_26_polcal1-22-03-11-13-25-59-0000.sci2', 'pixel_26_polcal1-22-03-11-13-25-59-0001.sci2'],'prefix':'2022-03-11_polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'} # BAD MEASUREMENT

# 2022-01-12
# {'pixel':5,'channel':25,'files':['pix_5_polcal1-22-01-12-12-17-03-0000.sci2','pix_5_polcal1-22-01-12-12-17-03-0001.sci2'],'prefix':'2022-01-12_polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':5,'channel':25,'files':['pix_5_polcal2-22-01-12-13-40-02-0000.sci2', 'pix_5_polcal2-22-01-12-13-40-02-0001.sci2'],'prefix':'2022-01-12_polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':19,'channel':23,'files':['pix_19_polcal1-22-01-12-11-32-02-0000.sci2', 'pix_19_polcal1-22-01-12-11-32-02-0001.sci2'],'prefix':'2022-01-12_polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
## {'pixel':19,'channel':23,'files':['pix_19_polcal2-22-01-12-13-26-00-0000.sci2', 'pix_19_polcal2-22-01-12-13-26-00-0001.sci2'],'prefix':'2022-01-12_polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'},\ ## BAD MEASUREMENT?
# {'pixel':23,'channel':24,'files':['pix_23_polcal1-22-01-12-11-47-09-0000.sci2', 'pix_23_polcal1-22-01-12-11-47-09-0001.sci2'],'prefix':'2022-01-12_polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':23,'channel':24,'files':['pix_23_polcal2-22-01-12-12-57-03-0000.sci2', 'pix_23_polcal2-22-01-12-12-57-03-0001.sci2'],'prefix':'2022-01-12_polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':26,'channel':21,'files':['pix_26_polcal1-22-01-12-12-03-03-0000.sci2', 'pix_26_polcal1-22-01-12-12-03-03-0001.sci2'],'prefix':'2022-01-12_polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':26,'channel':21,'files':['pix_26_polcal2-22-01-12-13-12-02-0000.sci2', 'pix_26_polcal2-22-01-12-13-12-02-0001.sci2'],'prefix':'2022-01-12_polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'}\

## 2021-11-23
{'pixel':41,'channel':4,'files':['pix41_polcal-21-11-22-12-25-54-0000.sci2','pix41_polcal-21-11-22-12-25-54-0001.sci2'],'prefix':'2021-11-21_polcal','fix_neg':False,'indir':basedir+'data/2021-11/'},\
{'pixel':41,'channel':4,'files':['PIX_41_polcal-21-11-23-12-23-03-0000.sci2','PIX_41_polcal-21-11-23-12-23-03-0001.sci2'],'prefix':'2021-11-21_polcal_2','fix_neg':False,'indir':basedir+'data/2021-11/'},\
{'pixel':42,'channel':6,'files':['PIX_42_polcal-21-11-23-12-08-09-0000.sci2','PIX_42_polcal-21-11-23-12-08-09-0001.sci2'],'prefix':'2021-11-21_polcal','fix_neg':False,'indir':basedir+'data/2021-11/'},\
{'pixel':63,'channel':5,'files':['PIX_63_polcal-21-11-23-11-47-05-0000.sci2','PIX_63_polcal-21-11-23-11-47-05-0001.sci2'],'prefix':'2021-11-21_polcal','fix_neg':False,'indir':basedir+'data/2021-11/'}


## Empty example
#{'pixel':,'channel':,'files'=,prefix=,fix_neg=False,'indir':basedir+'data/2022-01/'},\
]
for measurement in measurements:
	print('Pixel '+str(measurement['pixel']))
	print(measurement)
	try:
		data = run.get_sci(measurement['files'],quiet=True,indir=measurement['indir'])
		run.plot_sci(data,channels=[measurement['channel']],pixels=[measurement['pixel']],fix_neg=measurement['fix_neg'],offset=True,prefix=measurement['prefix'])
	except:
		input('Something went wrong, press return to continue')


# Run through the engineering data - note only the first file gets used anyway!
measurements = [\

## 2022-04-07
# {'pixel':42,'files':['pixel42_polcal1-22-04-07-14-24-44-0000.eng2','pixel42_polcal1-22-04-07-14-24-44-0001.eng2'],'prefix':'2022-04-07_polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':63,'files':['pixel63_polcal1-22-04-07-14-08-58-0000.eng2','pixel63_polcal1-22-04-07-14-08-58-0001.eng2'],'prefix':'2022-04-07_polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':41,'files':['pixel41_polcal1-22-04-07-13-54-58-0000.eng2','pixel41_polcal1-22-04-07-13-54-58-0001.eng2'],'prefix':'2022-04-07_polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal2-22-04-07-12-49-09-0000.eng2','pixel5_polcal2-22-04-07-12-49-09-0001.eng2'],'prefix':'2022-04-07_polcal2_new','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal2-22-04-07-11-57-01-0000.eng2','pixel5_polcal2-22-04-07-11-57-01-0001.eng2'],'prefix':'2022-04-07_polcal2','indir':basedir+'data/2022-04/'},\
# # NOTE: NEXT FILES MISNAMED, SAYS 26 ACTUALLY 19
# {'pixel':19,'files':['pixel26_polcal2-22-04-07-11-38-00-0000.eng2','pixel26_polcal2-22-04-07-11-38-00-0001.eng2'],'prefix':'2022-04-07_polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal2-22-04-07-11-21-57-0000.eng2','pixel26_polcal2-22-04-07-11-21-57-0001.eng2'],'prefix':'2022-04-07_polcal2_new','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal2-22-04-07-10-38-58-0000.eng2','pixel26_polcal2-22-04-07-10-38-58-0001.eng2'],'prefix':'2022-04-07_polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':23,'files':['pixel23_polcal2-22-04-07-09-46-00-0000.eng2','pixel23_polcal2-22-04-07-09-46-00-0001.eng2'],'prefix':'2022-04-07_polcal2','indir':basedir+'data/2022-04/'},\
# # {'pixel':23,'files':['pixel23_polcal2-22-04-07-09-43-59-0000.eng2'],'prefix':'2022-04-07_polcal2_bad','indir':basedir+'data/2022-04/'},\ # BAD MEASUREMENT

## 2022-04-06
# {'pixel':5,'files':['pixel5_polcal2-22-04-06-15-57-01-0000.eng2','pixel5_polcal2-22-04-06-15-57-01-0001.eng2'],'prefix':'2022-04-06_polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':19,'files':['pixel19_polcal2-22-04-06-15-42-59-0000.eng2','pixel19_polcal2-22-04-06-15-42-59-0001.eng2'],'prefix':'2022-04-06_polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal2-22-04-06-15-57-01-0000.eng2','pixel5_polcal2-22-04-06-15-57-01-0001.eng2'],'prefix':'2022-04-06_polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal2-22-04-06-15-19-58-0000.eng2','pixel26_polcal2-22-04-06-15-19-58-0001.eng2'],'prefix':'2022-04-06_polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':23,'files':['pixel23_polcal2-22-04-06-14-53-59-0000.eng2','pixel23_polcal2-22-04-06-14-53-59-0001.eng2'],'prefix':'2022-04-06_polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal1-22-04-06-13-55-56-0000.eng2','pixel5_polcal1-22-04-06-13-55-56-0001.eng2'],'prefix':'2022-04-06_polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':19,'files':['pixel19_polcal1-22-04-06-13-37-59-0000.eng2','pixel19_polcal1-22-04-06-13-37-59-0001.eng2'],'prefix':'2022-04-06_polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':23,'files':['pixel23_polcal1-22-04-06-11-57-59-0000.eng2','pixel23_polcal1-22-04-06-11-57-59-0001.eng2'],'prefix':'2022-04-06_polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal1-22-04-06-13-23-04-0000.eng2','pixel26_polcal1-22-04-06-13-23-04-0001.eng2'],'prefix':'2022-04-06_polcal1','indir':basedir+'data/2022-04/'}

## 2022-03-11
# {'pixel':5,'files':['pixel_5_polcal1-22-03-11-13-51-58-0000.eng2','pixel_5_polcal1-22-03-11-13-51-58-0001.eng2'],'prefix':'2022-03-11_polcal1','indir':basedir+'data/2022-03/'},\
# {'pixel':19,'files':['pixel_19_polcal1-22-03-11-13-38-00-0000.eng2','pixel_19_polcal1-22-03-11-13-38-00-0001.eng2'],'prefix':'2022-03-11_polcal1','indir':basedir+'data/2022-03/'},\
# {'pixel':23,'files':['pixel_23_polcal1-22-03-11-12-38-00-0000.eng2','pixel_23_polcal1-22-03-11-12-38-00-0001.eng2'],'prefix':'2022-03-11_polcal1_bad','indir':basedir+'data/2022-03/'},\
# {'pixel':23,'files':['pixel_23_polcal1-22-03-11-13-08-05-0000.eng2','pixel_23_polcal1-22-03-11-13-08-05-0001.eng2'],'prefix':'2022-03-11_polcal1','indir':basedir+'data/2022-03/'},\
# {'pixel':26,'files':['pixel_26_polcal1-22-03-11-13-23-06-0000.eng2','pixel_26_polcal1-22-03-11-13-23-06-0001.eng2'],'prefix':'2022-03-11_polcal1','indir':basedir+'data/2022-03/'}

## 2022-01-12
# {'pixel':5,'files':['pix_5_polcal1-22-01-12-12-14-04-0000.eng2','pix_5_polcal1-22-01-12-12-14-04-0001.eng2'],'prefix':'2022-01-12_polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':5,'files':['pix_5_polcal2-22-01-12-13-37-00-0000.eng2','pix_5_polcal2-22-01-12-13-37-00-0001.eng2'],'prefix':'2022-01-12_polcal2','indir':basedir+'data/2022-01/'},\
# {'pixel':19,'files':['pix_19_polcal1-22-01-12-11-28-06-0000.eng2','pix_19_polcal1-22-01-12-11-28-06-0001.eng2'],'prefix':'2022-01-12_polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':19,'files':['pix_19_polcal2-22-01-12-13-23-01-0000.eng2','pix_19_polcal2-22-01-12-13-23-01-0001.eng2'],'prefix':'2022-01-12_polcal2','indir':basedir+'data/2022-01/'},\
# {'pixel':23,'files':['pix_23_polcal1-22-01-12-11-44-03-0000.eng2','pix_23_polcal1-22-01-12-11-44-03-0001.eng2'],'prefix':'2022-01-12_polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':23,'files':['pix_23_polcal2-22-01-12-12-54-02-0000.eng2','pix_23_polcal2-22-01-12-12-54-02-0001.eng2'],'prefix':'2022-01-12_polcal2','indir':basedir+'data/2022-01/'},\
# {'pixel':26,'files':['pix_26_polcal1-22-01-12-12-00-09-0000.eng2','pix_26_polcal1-22-01-12-12-00-09-0001.eng2'],'prefix':'2022-01-12_polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':26,'files':['pix_26_polcal2-22-01-12-13-09-04-0000.eng2','pix_26_polcal2-22-01-12-13-09-04-0001.eng2'],'prefix':'2022-01-12_polcal2','indir':basedir+'data/2022-01/'}\

## 2021-11-22
{'pixel':41,'files':['pix41_polcal-21-11-22-12-21-19-0000.eng2','pix41_polcal-21-11-22-12-21-19-0001.eng2'],'prefix':'2021-11-22_polcal','indir':basedir+'data/2021-11/'},\
{'pixel':41,'files':['PIX_41_polcal-21-11-23-12-20-01-0000.eng2','PIX_41_polcal-21-11-23-12-20-01-0001.eng2'],'prefix':'2021-11-23_polcal_2','indir':basedir+'data/2021-11/'},\
{'pixel':63,'files':['PIX_63_polcal-21-11-23-11-44-06-0000.eng2','PIX_63_polcal-21-11-23-11-44-06-0001.eng2'],'prefix':'2021-11-23','indir':basedir+'data/2021-11/'},\
{'pixel':42,'files':['PIX_42_polcal-21-11-23-12-04-01-0000.eng2','PIX_42_polcal-21-11-23-12-04-01-0001.eng2'],'prefix':'2021-11-23','indir':basedir+'data/2021-11/'},\
{'pixel':23,'files':['PIX_23_polcal-21-11-23-13-32-01-0000.eng2','PIX_23_polcal-21-11-23-13-32-01-0001.eng2'],'prefix':'2021-11-23','indir':basedir+'data/2021-11/'},\

## Empty example
# {'pixel':,'files'=,prefix=,'indir':basedir+'data/2022-01/'},\
]
for measurement in measurements:
	print('Pixel '+str(measurement['pixel']))
	print(measurement)
	try:
		data = run.get_eng(measurement['files'],quiet=True,indir=measurement['indir'])
		run.plot_eng(data,pixel=measurement['pixel'],prefix=measurement['prefix'])
	except:
		input('Something went wrong, press return to continue')
