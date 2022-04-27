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
basedir = '/Volumes/WD12TB/quijote2/'
# basedir = '/Users/mpeel/Documents/git/tfgi/'
# basedir = '/Volumes/proyectos-1/quijote2/'

# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
# outdir = '/Volumes/WD12TB/quijote2/output/'
outdir = '/Users/mpeel/Documents/git/tfgi/'

# Start the class
# NOTE: this code does not actually use the pixel masterfiles, since things might change during a measurement. They are hard-coded for the sci data (they are not needed for eng)
run = tfgi.tfgi(outdir=outdir,\
	datadir=basedir+'tod/',\
	pixelfileloc=basedir+'etc/qt2_pixel_masterfile.',\
	pixelposfileloc=basedir+'etc/tgi_fgi_horn_positions_table.txt',\
	polcalfileloc=basedir+'etc/qt2_masterfiles_new/qt2_polcal.')

# Still data to add to this list from 2018-06 and 2018-05 and 2017-03.data and 2017-03 and 2017-02 and 2017-01 and 2016-12 and 2016-11 and 2016-10

# Run through the scientific data
# NOTE: if these don't run, make sure you have a ,\ at the end of each line - EXCEPT the last one, where it should end with a }
measurements = [\

## 2022-04-08
# {'pixel':41,'channel':4,'files':['pixel41_polcal1-22-04-08-10-16-58-0000.sci2','pixel41_polcal1-22-04-08-10-16-58-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':63,'channel':5,'files':['pixel63_polcal1-22-04-08-10-31-59-0000.sci2','pixel63_polcal1-22-04-08-10-31-59-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':42,'channel':6,'files':['pixel42_polcal1-22-04-08-10-47-02-0000.sci2','pixel42_polcal1-22-04-08-10-47-02-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':41,'channel':4,'files':['pixel41_polcal1-22-04-08-11-28-00-0000.sci2','pixel41_polcal1-22-04-08-11-28-00-0001.sci2'],'prefix':'polcal1_repeat','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':63,'channel':5,'files':['pixel63_polcal1-22-04-08-11-50-00-0000.sci2','pixel63_polcal1-22-04-08-11-50-00-0001.sci2'],'prefix':'polcal1_repeat','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':42,'channel':6,'files':['pixel42_polcal1-22-04-08-12-06-58-0000.sci2','pixel42_polcal1-22-04-08-12-06-58-0001.sci2'],'prefix':'polcal1_repeat2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':41,'channel':4,'files':['pixel41_polcal1-22-04-08-12-36-00-0000.sci2','pixel41_polcal1-22-04-08-12-36-00-0001.sci2'],'prefix':'polcal1_repeat2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':41,'channel':4,'files':['pixel41_polcal1-22-04-08-13-00-58-0000.sci2','pixel41_polcal1-22-04-08-13-00-58-0001.sci2'],'prefix':'polcal1_repeat3','fix_neg':False,'indir':basedir+'data/2022-04/'},\

## 2022-04-07
# {'pixel':42,'channel':6,'files':['pixel42_polcal1-22-04-07-14-28-23-0000.sci2','pixel42_polcal1-22-04-07-14-28-23-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':63,'channel':5,'files':['pixel63_polcal1-22-04-07-14-12-03-0000.sci2','pixel63_polcal1-22-04-07-14-12-03-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':41,'channel':4,'files':['pixel41_polcal1-22-04-07-13-57-58-0000.sci2','pixel41_polcal1-22-04-07-13-57-58-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':5,'channel':25,'files':['pixel5_polcal2-22-04-07-12-52-03-0000.sci2','pixel5_polcal2-22-04-07-12-52-03-0001.sci2'],'prefix':'polcal2_new','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':5,'channel':25,'files':['pixel5_polcal2-22-04-07-11-59-59-0000.sci2','pixel5_polcal2-22-04-07-11-59-59-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # NOTE: NEXT FILES MISNAMED, SAYS 26 ACTUALLY 19
# {'pixel':19,'channel':23,'files':['pixel26_polcal2-22-04-07-11-40-59-0000.sci2','pixel26_polcal2-22-04-07-11-40-59-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # LAST BAD MEASUREMENT?
# {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-07-11-24-59-0000.sci2','pixel26_polcal2-22-04-07-11-24-59-0001.sci2'],'prefix':'polcal2_new','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-07-10-42-59-0000.sci2','pixel26_polcal2-22-04-07-10-42-59-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # LAST BAD MEASUREMENT?
# {'pixel':23,'channel':24,'files':['pixel23_polcal2-22-04-07-09-48-58-0000.sci2','pixel23_polcal2-22-04-07-09-48-58-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\

## 2022-04-06
# {'pixel':5,'channel':25,'files':['pixel5_polcal2-22-04-06-16-03-58-0000.sci2','pixel5_polcal2-22-04-06-16-03-58-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':19,'channel':23,'files':['pixel19_polcal2-22-04-06-15-46-00-0000.sci2','pixel19_polcal2-22-04-06-15-46-00-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # LAST MEASUREMENT BAD???
# {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-06-15-31-58-0000.sci2','pixel26_polcal2-22-04-06-15-31-58-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # LAST MEASUREMENT BAD??
# {'pixel':26,'channel':21,'files':['pixel26_polcal2-22-04-06-15-24-02-0000.sci2','pixel26_polcal2-22-04-06-15-24-02-0001.sci2'],'prefix':'polcal2_bad','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# # LAST MEASUREMENT BAD??
# {'pixel':23,'channel':24,'files':['pixel23_polcal2-22-04-06-15-06-58-0000.sci2','pixel23_polcal2-22-04-06-15-06-58-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':23,'channel':24,'files':['pixel23_polcal2-22-04-06-14-58-58-0000.sci2','pixel23_polcal2-22-04-06-14-58-58-0001.sci2'],'prefix':'polcal2_bad','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':5,'channel':25,'files':['pixel5_polcal1-22-04-06-13-59-57-0000.sci2','pixel5_polcal1-22-04-06-13-59-57-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':19,'channel':23,'files':['pixel19_polcal1-22-04-06-13-40-56-0000.sci2','pixel19_polcal1-22-04-06-13-40-56-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':26,'channel':21,'files':['pixel26_polcal1-22-04-06-13-26-00-0000.sci2','pixel26_polcal1-22-04-06-13-26-00-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\
# {'pixel':23,'channel':24,'files':['pixel23_polcal1-22-04-06-12-08-00-0000.sci2','pixel23_polcal1-22-04-06-12-08-00-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-04/'},\

## 2022-03-11
# {'pixel':5,'channel':25,'files':['pixel_5_polcal1-22-03-11-13-55-19-0000.sci2','pixel_5_polcal1-22-03-11-13-55-19-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-03/'},\
# {'pixel':19,'channel':23,'files':['pixel_19_polcal1-22-03-11-13-41-11-0000.sci2', 'pixel_19_polcal1-22-03-11-13-41-11-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-03/'},\
# {'pixel':23,'channel':24,'files':['pixel_23_polcal1-22-03-11-13-11-08-0000.sci2', 'pixel_23_polcal1-22-03-11-13-11-08-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-03/'},\
# {'pixel':26,'channel':21,'files':['pixel_26_polcal1-22-03-11-13-25-59-0000.sci2', 'pixel_26_polcal1-22-03-11-13-25-59-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# # LAST BAD MEASUREMENT

# 2022-01-12
# {'pixel':5,'channel':25,'files':['pix_5_polcal1-22-01-12-12-17-03-0000.sci2','pix_5_polcal1-22-01-12-12-17-03-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':5,'channel':25,'files':['pix_5_polcal2-22-01-12-13-40-02-0000.sci2', 'pix_5_polcal2-22-01-12-13-40-02-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':19,'channel':23,'files':['pix_19_polcal1-22-01-12-11-32-02-0000.sci2', 'pix_19_polcal1-22-01-12-11-32-02-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':19,'channel':23,'files':['pix_19_polcal2-22-01-12-13-26-00-0000.sci2', 'pix_19_polcal2-22-01-12-13-26-00-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# ## LAST BAD MEASUREMENT?
# {'pixel':23,'channel':24,'files':['pix_23_polcal1-22-01-12-11-47-09-0000.sci2', 'pix_23_polcal1-22-01-12-11-47-09-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':23,'channel':24,'files':['pix_23_polcal2-22-01-12-12-57-03-0000.sci2', 'pix_23_polcal2-22-01-12-12-57-03-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':26,'channel':21,'files':['pix_26_polcal1-22-01-12-12-03-03-0000.sci2', 'pix_26_polcal1-22-01-12-12-03-03-0001.sci2'],'prefix':'polcal1','fix_neg':False,'indir':basedir+'data/2022-01/'},\
# {'pixel':26,'channel':21,'files':['pix_26_polcal2-22-01-12-13-12-02-0000.sci2', 'pix_26_polcal2-22-01-12-13-12-02-0001.sci2'],'prefix':'polcal2','fix_neg':False,'indir':basedir+'data/2022-01/'},\

## 2021-11-23
# {'pixel':41,'channel':4,'files':['pix41_polcal-21-11-22-12-25-54-0000.sci2','pix41_polcal-21-11-22-12-25-54-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2021-11/'},\
# {'pixel':41,'channel':4,'files':['PIX_41_polcal-21-11-23-12-23-03-0000.sci2','PIX_41_polcal-21-11-23-12-23-03-0001.sci2'],'prefix':'polcal_2','fix_neg':False,'indir':basedir+'data/2021-11/'},\
# {'pixel':42,'channel':6,'files':['PIX_42_polcal-21-11-23-12-08-09-0000.sci2','PIX_42_polcal-21-11-23-12-08-09-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2021-11/'},\
# {'pixel':63,'channel':5,'files':['PIX_63_polcal-21-11-23-11-47-05-0000.sci2','PIX_63_polcal-21-11-23-11-47-05-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2021-11/'},\

## 2019-04-10
# {'pixel':41,'channel':4,'files':['Pix41_cal-19-04-10-11-31-14-0000.sci2','Pix41_cal-19-04-10-11-31-14-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# {'pixel':42,'channel':6,'files':['Pix42_cal-19-04-10-12-07-14-0000.sci2','Pix42_cal-19-04-10-12-07-14-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# {'pixel':63,'channel':5,'files':['Pix63_cal-19-04-10-11-48-14-0000.sci2','Pix63_cal-19-04-10-11-48-14-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2019-04/'},\

## 2019-04-09
# {'pixel':5,'channel':25,'files':['pix5_cal-19-04-09-13-14-17-0000.sci2', 'pix5_cal-19-04-09-13-14-17-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# {'pixel':17,'channel':23,'files':['pix17_cal-19-04-09-12-56-16-0000.sci2', 'pix17_cal-19-04-09-12-56-16-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# {'pixel':23,'channel':24,'files':['pix23_cal-19-04-09-13-33-21-0000.sci2', 'pix23_cal-19-04-09-13-33-21-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# {'pixel':26,'channel':21,'files':['pix26_cal-19-04-09-11-49-35-0000.sci2', 'pix26_cal-19-04-09-11-49-35-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# # LAST BAD MEASUREMENT?
# {'pixel':26,'channel':21,'files':['pix26_cal-19-04-09-12-20-15-0000.sci2', 'pix26_cal-19-04-09-12-20-15-0001.sci2'],'prefix':'polcal_2','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# # LAST BAD MEASUREMENT?
# {'pixel':26,'channel':21,'files':['pix26_cal-19-04-09-12-36-17-0000.sci2', 'pix26_cal-19-04-09-12-36-17-0001.sci2'],'prefix':'polcal_3','fix_neg':False,'indir':basedir+'data/2019-04/'},\
# # LAST BAD MEASUREMENT?

## 2018-09
# {'pixel':2,'channel':20,'files':['Pix_2-18-09-12-11-44-10-0000.sci2','Pix_2-18-09-12-11-44-10-0001.sci2', 'Pix_2-18-09-12-11-44-10-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':5,'channel':25,'files':['Pix_5-18-09-18-10-57-45-0000.sci2','Pix_5-18-09-18-10-57-45-0001.sci2', 'Pix_5-18-09-18-10-57-45-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':8,'channel':26,'files':['Pix_8-18-09-18-12-42-37-0000.sci2','Pix_8-18-09-18-12-42-37-0001.sci2', 'Pix_8-18-09-18-12-42-37-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':10,'channel':30,'files':['Pix_10-18-09-18-12-20-38-0000.sci2','Pix_10-18-09-18-12-20-38-0001.sci2', 'Pix_10-18-09-18-12-20-38-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':12,'channel':28,'files':['Pix_12-18-09-18-11-53-36-0000.sci2','Pix_12-18-09-18-11-53-36-0001.sci2', 'Pix_12-18-09-18-11-53-36-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':15,'channel':29,'files':['Pix_15-18-09-18-10-30-37-0000.sci2','Pix_15-18-09-18-10-30-37-0001.sci2', 'Pix_15-18-09-18-10-30-37-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':17,'channel':23,'files':['Pix_17-18-09-18-09-52-36-0000.sci2','Pix_17-18-09-18-09-52-36-0001.sci2', 'Pix_17-18-09-18-09-52-36-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':19,'channel':22,'files':['Pix_19-18-09-12-12-12-36-0000.sci2','Pix_19-18-09-12-12-12-36-0001.sci2', 'Pix_19-18-09-12-12-12-36-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':20,'channel':27,'files':['Pix_20-18-09-18-11-25-37-0000.sci2','Pix_20-18-09-18-11-25-37-0001.sci2', 'Pix_20-18-09-18-11-25-37-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':23,'channel':24,'files':['Pix_23-18-09-18-10-10-37-0000.sci2','Pix_23-18-09-18-10-10-37-0001.sci2', 'Pix_23-18-09-18-10-10-37-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':27,'channel':23,'files':['Pix_27-18-09-11-12-44-36-0000.sci2','Pix_27-18-09-11-12-44-36-0001.sci2', 'Pix_27-18-09-11-12-44-36-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':42,'channel':6,'files':['Pix_42-18-09-11-12-25-34-0000.sci2','Pix_42-18-09-11-12-25-34-0001.sci2', 'Pix_42-18-09-11-12-25-34-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':48,'channel':13,'files':['Pix_48-18-09-11-10-29-44-0000.sci2','Pix_48-18-09-11-10-29-44-0001.sci2', 'Pix_48-18-09-11-10-29-44-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':50,'channel':7,'files':['Pix_50-18-09-11-13-48-35-0000.sci2','Pix_50-18-09-11-13-48-35-0001.sci2', 'Pix_50-18-09-11-13-48-35-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':51,'channel':9,'files':['Pix_51-18-09-11-12-09-37-0000.sci2','Pix_51-18-09-11-12-09-37-0001.sci2', 'Pix_51-18-09-11-12-09-37-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':53,'channel':14,'files':['Pix_53-18-09-11-10-04-36-0000.sci2','Pix_53-18-09-11-10-04-36-0001.sci2', 'Pix_53-18-09-11-10-04-36-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':55,'channel':10,'files':['Pix_55-18-09-11-14-28-37-0000.sci2','Pix_55-18-09-11-14-28-37-0001.sci2', 'Pix_55-18-09-11-14-28-37-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':63,'channel':5,'files':['Pix_63-18-09-11-10-46-37-0000.sci2','Pix_63-18-09-11-10-46-37-0001.sci2', 'Pix_63-18-09-11-10-46-37-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':64,'channel':12,'files':['Pix_64-18-09-11-11-45-54-0000.sci2','Pix_64-18-09-11-11-45-54-0001.sci2', 'Pix_64-18-09-11-11-45-54-0002.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':65,'channel':8,'files':['Pix_65-18-09-11-13-16-37-0000.sci2','Pix_65-18-09-11-13-16-37-0001.sci2'],'prefix':'polcal','fix_neg':False,'indir':basedir+'data/2018-09/'}

## 2018-07  NOT YET RUN
# {'pixel':44,'channel':0,'files':['Pix44_cal-18-07-02-10-56-29-0000.sci2','Pix44_cal-18-07-02-10-56-29-0001.sci2', 'Pix44_cal-18-07-02-10-56-29-0002.sci2','Pix44_cal-18-07-02-10-56-29-0003.sci2','Pix44_cal-18-07-02-10-56-29-0004.sci2','Pix44_cal-18-07-02-10-56-29-0005.sci2'],'prefix':'2018-07_polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':64,'channel':0,'files':['Pix64_cal-18-07-02-12-38-28-0000.sci2','Pix64_cal-18-07-02-12-38-28-0001.sci2', 'Pix64_cal-18-07-02-12-38-28-0002.sci2','Pix64_cal-18-07-02-12-38-28-0003.sci2','Pix64_cal-18-07-02-12-38-28-0004.sci2','Pix64_cal-18-07-02-12-38-28-0005.sci2'],'prefix':'2018-07_polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\
# {'pixel':65,'channel':0,'files':['Pix65_cal-18-07-02-13-54-27-0000.sci2','Pix65_cal-18-07-02-13-54-27-0001.sci2', 'Pix65_cal-18-07-02-13-54-27-0002.sci2','Pix65_cal-18-07-02-13-54-27-0003.sci2','Pix65_cal-18-07-02-13-54-27-0004.sci2','Pix65_cal-18-07-02-13-54-27-0005.sci2'],'prefix':'2018-07_polcal','fix_neg':False,'indir':basedir+'data/2018-09/'},\


## Empty example
#{'pixel':,'channel':,'files'=,prefix=,fix_neg=False,'indir':basedir+'data/2022-01/'},\
]
for measurement in measurements:
	print('Pixel '+str(measurement['pixel']))
	print(measurement)
	try:
		data = run.get_sci(measurement['files'],quiet=True,indir=measurement['indir'])
		temp = measurement['files'][0].split('-')
		date = '20' + temp[1] + '-' + temp[2] + '-' + temp[3] + '_' + temp[4] + '-' + temp[5]
		run.plot_sci(data,channels=[measurement['channel']],pixels=[measurement['pixel']],fix_neg=measurement['fix_neg'],offset=True,prefix=date+'_'+measurement['prefix'])
	except:
		input('Something went wrong, press return to continue')

# Run through the engineering data - note only the first file gets used anyway!
measurements = [\

## 2022-04-08
# {'pixel':41,'files':['pixel41_polcal1-22-04-08-10-14-00-0000.eng2','pixel41_polcal1-22-04-08-10-14-00-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':63,'files':['pixel63_polcal1-22-04-08-10-29-01-0000.eng2','pixel63_polcal1-22-04-08-10-29-01-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':42,'files':['pixel42_polcal1-22-04-08-10-44-02-0001.eng2','pixel42_polcal1-22-04-08-10-44-02-0000.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':41,'files':['pixel41_polcal1-22-04-08-11-19-59-0000.eng2','pixel41_polcal1-22-04-08-11-19-59-0001.eng2'],'prefix':'polcal1_repeat','indir':basedir+'data/2022-04/'},\
# {'pixel':63,'files':['pixel63_polcal1-22-04-08-11-38-58-0000.eng2','pixel63_polcal1-22-04-08-11-38-58-0001.eng2'],'prefix':'polcal1_repeat','indir':basedir+'data/2022-04/'},\
# {'pixel':42,'files':['pixel42_polcal1-22-04-08-10-44-02-0001.eng2','pixel42_polcal1-22-04-08-10-44-02-0000.eng2'],'prefix':'polcal1_repeat','indir':basedir+'data/2022-04/'},\
# {'pixel':42,'files':['pixel42_polcal1-22-04-08-12-04-01-0000.eng2','pixel42_polcal1-22-04-08-12-04-01-0001.eng2'],'prefix':'polcal1_repeat2','indir':basedir+'data/2022-04/'},\
# {'pixel':41,'files':['pixel41_polcal1-22-04-08-12-33-03-0000.eng2','pixel41_polcal1-22-04-08-12-33-03-0001.eng2'],'prefix':'polcal1_repeat2','indir':basedir+'data/2022-04/'},\
# {'pixel':41,'files':['pixel41_polcal1-22-04-08-12-57-59-0000.eng2','pixel41_polcal1-22-04-08-12-57-59-0001.eng2'],'prefix':'polcal1_repeat3','indir':basedir+'data/2022-04/'},\
#
# ## 2022-04-07
# {'pixel':42,'files':['pixel42_polcal1-22-04-07-14-24-44-0000.eng2','pixel42_polcal1-22-04-07-14-24-44-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':63,'files':['pixel63_polcal1-22-04-07-14-08-58-0000.eng2','pixel63_polcal1-22-04-07-14-08-58-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':41,'files':['pixel41_polcal1-22-04-07-13-54-58-0000.eng2','pixel41_polcal1-22-04-07-13-54-58-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal2-22-04-07-12-49-09-0000.eng2','pixel5_polcal2-22-04-07-12-49-09-0001.eng2'],'prefix':'polcal2_new','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal2-22-04-07-11-57-01-0000.eng2','pixel5_polcal2-22-04-07-11-57-01-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# # NOTE: NEXT FILES MISNAMED, SAYS 26 ACTUALLY 19
# {'pixel':19,'files':['pixel26_polcal2-22-04-07-11-38-00-0000.eng2','pixel26_polcal2-22-04-07-11-38-00-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal2-22-04-07-11-21-57-0000.eng2','pixel26_polcal2-22-04-07-11-21-57-0001.eng2'],'prefix':'polcal2_new','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal2-22-04-07-10-38-58-0000.eng2','pixel26_polcal2-22-04-07-10-38-58-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':23,'files':['pixel23_polcal2-22-04-07-09-46-00-0000.eng2','pixel23_polcal2-22-04-07-09-46-00-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':23,'files':['pixel23_polcal2-22-04-07-09-43-59-0000.eng2'],'prefix':'polcal2_bad','indir':basedir+'data/2022-04/'},\
# # LAST BAD MEASUREMENT
#
# ## 2022-04-06
# {'pixel':5,'files':['pixel5_polcal2-22-04-06-15-57-01-0000.eng2','pixel5_polcal2-22-04-06-15-57-01-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':19,'files':['pixel19_polcal2-22-04-06-15-42-59-0000.eng2','pixel19_polcal2-22-04-06-15-42-59-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal2-22-04-06-15-57-01-0000.eng2','pixel5_polcal2-22-04-06-15-57-01-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal2-22-04-06-15-19-58-0000.eng2','pixel26_polcal2-22-04-06-15-19-58-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':23,'files':['pixel23_polcal2-22-04-06-14-53-59-0000.eng2','pixel23_polcal2-22-04-06-14-53-59-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-04/'},\
# {'pixel':5,'files':['pixel5_polcal1-22-04-06-13-55-56-0000.eng2','pixel5_polcal1-22-04-06-13-55-56-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':19,'files':['pixel19_polcal1-22-04-06-13-37-59-0000.eng2','pixel19_polcal1-22-04-06-13-37-59-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':23,'files':['pixel23_polcal1-22-04-06-11-57-59-0000.eng2','pixel23_polcal1-22-04-06-11-57-59-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
# {'pixel':26,'files':['pixel26_polcal1-22-04-06-13-23-04-0000.eng2','pixel26_polcal1-22-04-06-13-23-04-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-04/'},\
#
# ## 2022-03-11
# {'pixel':5,'files':['pixel_5_polcal1-22-03-11-13-51-58-0000.eng2','pixel_5_polcal1-22-03-11-13-51-58-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-03/'},\
# {'pixel':19,'files':['pixel_19_polcal1-22-03-11-13-38-00-0000.eng2','pixel_19_polcal1-22-03-11-13-38-00-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-03/'},\
# {'pixel':23,'files':['pixel_23_polcal1-22-03-11-12-38-00-0000.eng2','pixel_23_polcal1-22-03-11-12-38-00-0001.eng2'],'prefix':'polcal1_bad','indir':basedir+'data/2022-03/'},\
# {'pixel':23,'files':['pixel_23_polcal1-22-03-11-13-08-05-0000.eng2','pixel_23_polcal1-22-03-11-13-08-05-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-03/'},\
# {'pixel':26,'files':['pixel_26_polcal1-22-03-11-13-23-06-0000.eng2','pixel_26_polcal1-22-03-11-13-23-06-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-03/'},\
#
# ## 2022-01-12
# {'pixel':5,'files':['pix_5_polcal1-22-01-12-12-14-04-0000.eng2','pix_5_polcal1-22-01-12-12-14-04-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':5,'files':['pix_5_polcal2-22-01-12-13-37-00-0000.eng2','pix_5_polcal2-22-01-12-13-37-00-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-01/'},\
# {'pixel':19,'files':['pix_19_polcal1-22-01-12-11-28-06-0000.eng2','pix_19_polcal1-22-01-12-11-28-06-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':19,'files':['pix_19_polcal2-22-01-12-13-23-01-0000.eng2','pix_19_polcal2-22-01-12-13-23-01-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-01/'},\
# {'pixel':23,'files':['pix_23_polcal1-22-01-12-11-44-03-0000.eng2','pix_23_polcal1-22-01-12-11-44-03-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':23,'files':['pix_23_polcal2-22-01-12-12-54-02-0000.eng2','pix_23_polcal2-22-01-12-12-54-02-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-01/'},\
# {'pixel':26,'files':['pix_26_polcal1-22-01-12-12-00-09-0000.eng2','pix_26_polcal1-22-01-12-12-00-09-0001.eng2'],'prefix':'polcal1','indir':basedir+'data/2022-01/'},\
# {'pixel':26,'files':['pix_26_polcal2-22-01-12-13-09-04-0000.eng2','pix_26_polcal2-22-01-12-13-09-04-0001.eng2'],'prefix':'polcal2','indir':basedir+'data/2022-01/'},\

## 2021-11-22
{'pixel':41,'files':['pix41_polcal-21-11-22-12-21-19-0000.eng2','pix41_polcal-21-11-22-12-21-19-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2021-11/'},\
{'pixel':41,'files':['PIX_41_polcal-21-11-23-12-20-01-0000.eng2','PIX_41_polcal-21-11-23-12-20-01-0001.eng2'],'prefix':'polcal_2','indir':basedir+'data/2021-11/'},\
{'pixel':63,'files':['PIX_63_polcal-21-11-23-11-44-06-0000.eng2','PIX_63_polcal-21-11-23-11-44-06-0001.eng2'],'prefix':'','indir':basedir+'data/2021-11/'},\
{'pixel':42,'files':['PIX_42_polcal-21-11-23-12-04-01-0000.eng2','PIX_42_polcal-21-11-23-12-04-01-0001.eng2'],'prefix':'','indir':basedir+'data/2021-11/'},\
{'pixel':23,'files':['PIX_23_polcal-21-11-23-13-32-01-0000.eng2','PIX_23_polcal-21-11-23-13-32-01-0001.eng2'],'prefix':'','indir':basedir+'data/2021-11/'},\

## 2019-04-10
{'pixel':41,'files':['Pix41_cal-19-04-10-11-26-25-0000.eng2','Pix41_cal-19-04-10-11-26-25-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2019-04/'},\
{'pixel':42,'files':['Pix42_cal-19-04-10-12-03-14-0000.eng2','Pix42_cal-19-04-10-12-03-14-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2019-04/'},\
{'pixel':63,'files':['Pix63_cal-19-04-10-11-43-50-0000.eng2','Pix63_cal-19-04-10-11-43-50-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2019-04/'},\

## 2019-04-09
{'pixel':5,'files':['pix5_cal-19-04-09-13-10-07-0000.eng2','pix5_cal-19-04-09-13-10-07-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2019-04/'},\
{'pixel':17,'files':['pix17_cal-19-04-09-12-51-20-0000.eng2','pix17_cal-19-04-09-12-51-20-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2019-04/'},\
{'pixel':23,'files':['pix23_cal-19-04-09-13-28-32-0000.eng2','pix23_cal-19-04-09-13-28-32-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2019-04/'},\
{'pixel':26,'files':['pix26_cal-19-04-09-11-45-34-0000.eng2','pix26_cal-19-04-09-11-45-34-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2019-04/'},\
{'pixel':26,'files':['pix26_cal-19-04-09-11-48-15-0000.eng2'],'prefix':'2021-04-09_polcal_2','indir':basedir+'data/2019-04/'},\
# LAST BAD MEASUREMENT?

## 2018-09
{'pixel':2,'files':['Pix_2-18-09-12-11-41-37-0000.eng2','Pix_2-18-09-12-11-41-37-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':5,'files':['Pix_5-18-09-18-10-54-37-0000.eng2','Pix_5-18-09-18-10-54-37-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':8,'files':['Pix_8-18-09-18-12-39-38-0000.eng2','Pix_8-18-09-18-12-39-38-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':10,'files':['Pix_10-18-09-18-12-15-46-0000.eng2','Pix_10-18-09-18-12-15-46-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':12,'files':['Pix_12-18-09-18-11-50-40-0000.eng2','Pix_12-18-09-18-11-50-40-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':15,'files':['Pix_15-18-09-18-10-27-37-0000.eng2','Pix_15-18-09-18-10-27-37-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':17,'files':['Pix_17-18-09-18-09-49-38-0000.eng2','Pix_17-18-09-18-09-49-38-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':19,'files':['Pix_19-18-09-12-12-09-35-0000.eng2','Pix_19-18-09-12-12-09-35-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':20,'files':['Pix_20-18-09-18-11-22-38-0000.eng2','Pix_20-18-09-18-11-22-38-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':23,'files':['Pix_23-18-09-18-10-07-38-0000.eng2','Pix_23-18-09-18-10-07-38-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':27,'files':['Pix_27-18-09-11-12-41-48-0000.eng2','Pix_27-18-09-11-12-41-48-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':42,'files':['Pix_42-18-09-11-12-22-47-0000.eng2','Pix_42-18-09-11-12-22-47-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':48,'files':['Pix_48-18-09-11-10-26-45-0000.eng2','Pix_48-18-09-11-10-26-45-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':50,'files':['Pix_50-18-09-11-13-45-39-0000.eng2','Pix_50-18-09-11-13-45-39-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':51,'files':['Pix_51-18-09-11-12-06-35-0000.eng2','Pix_51-18-09-11-12-06-35-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':53,'files':['Pix_53-18-09-11-10-01-35-0000.eng2','Pix_53-18-09-11-10-01-35-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':55,'files':['Pix_55-18-09-11-14-25-34-0000.eng2','Pix_55-18-09-11-14-25-34-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':63,'files':['Pix_63-18-09-11-10-43-35-0000.eng2','Pix_63-18-09-11-10-43-35-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':64,'files':['Pix_64-18-09-11-11-42-12-0000.eng2','Pix_64-18-09-11-11-42-12-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\
{'pixel':65,'files':['Pix_65-18-09-11-13-13-48-0000.eng2','Pix_65-18-09-11-13-13-48-0001.eng2'],'prefix':'polcal','indir':basedir+'data/2018-09/'},\

## 2018-07
{'pixel':44,'files':['Pix44_cal-18-07-02-10-34-30-0000.eng2'],'prefix':'2018-07_polcal','indir':basedir+'data/2018-07/'},\
# LAST BAD MEASUREMENT?
{'pixel':44,'files':['Pix44_cal-18-07-02-10-53-28-0000.eng2','Pix44_cal-18-07-02-10-53-28-0001.eng2','Pix44_cal-18-07-02-10-53-28-0002.eng2'],'prefix':'polcal_2','indir':basedir+'data/2018-07/'},\
{'pixel':64,'files':['Pix64_cal-18-07-02-12-35-27-0000.eng2','Pix64_cal-18-07-02-12-35-27-0001.eng2','Pix64_cal-18-07-02-12-35-27-0002.eng2'],'prefix':'polcal','indir':basedir+'data/2018-07/'},\
{'pixel':65,'files':['Pix65_cal-18-07-02-13-51-29-0000.eng2','Pix65_cal-18-07-02-13-51-29-0001.eng2','Pix65_cal-18-07-02-13-51-29-0002.eng2'],'prefix':'polcal','indir':basedir+'data/2018-07/'},\

## Empty example
# {'pixel':,'files'=,prefix=,'indir':basedir+'data/2022-01/'},\
]
for measurement in measurements:
	print('Pixel '+str(measurement['pixel']))
	print(measurement)
	try:
		data = run.get_eng(measurement['files'],quiet=True,indir=measurement['indir'])
		temp = measurement['files'][0].split('-')
		date = '20' + temp[1] + '-' + temp[2] + '-' + temp[3] + '_' + temp[4] + '-' + temp[5]
		run.plot_eng(data,pixel=measurement['pixel'],prefix=date+'_'+measurement['prefix'])
	except:
		input('Something went wrong, press return to continue')
