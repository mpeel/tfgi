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

#
# NB: documentation is at:
# https://docs.google.com/document/d/13CzBq_ABo7NTrcoM9C_v25ws_2t3GHdxKEKEkGspegA
#

# This is needed if you want to write out lots of plots
mpl.rcParams['agg.path.chunksize'] = 10000
plt.ioff()

# Set this to where you have a copy of the data
basedir = '/Volumes/WD12TB/quijote2/'
# Set this to where you want to output stuff. The directory should already exist.
# NB: subdirectories will automatically be created for each dataset.
outdir = '/Volumes/WD12TB/quijote2/output/'
# Start the class
run = tfgi.tfgi(outdir=outdir,\
	datadir=basedir+'tod/',\
	pixelfileloc=basedir+'etc/qt2_pixel_masterfile.',\
	pixelposfileloc=basedir+'etc/tgi_fgi_horn_positions_table.txt',\
	polcalfileloc=basedir+'etc/qt2_polcal.',
	nside = 512)

# Search for CRAB and MOON observations in April 2019, and analyse them.
# datasets1 = run.find_observations('CASS-211115-1805')
# datasets2 = run.find_observations('CASS-211115-1909')
# datasets3 = run.find_observations('CASS-211115-2039')
# datasets4 = run.find_observations('CRAB-211115-2250')
datasets = run.find_observations('LOCALMAP-211115')
# datasets = run.find_observations('MOON-2111')
# # datasets2 = run.find_observations('DIP000')
# # datasets = list(set(datasets1) | set(datasets2))
# datasets = list(set(datasets1) | set(datasets2) | set(datasets3) | set(datasets4))
# datasets1 = run.find_observations('DIP000-211117')
# datasets2 = run.find_observations('CASS-211117')
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
	run.analyse_tod(dataset,plotlimit=0.001,dopol=True,plottods=True)

exit()
