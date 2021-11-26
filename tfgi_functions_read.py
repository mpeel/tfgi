#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Various functions related to reading tgi data, separate from the class
#
# Version history:
#
# 02-May-2019  M. Peel       Split from tgi.py

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import pandas as pd
import scipy.fftpack
from scipy import signal, optimize
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
from astropy.time import Time
from astropy.coordinates import Angle
from scipy.optimize import curve_fit
import os
import astroplan

# Read in a single pixel masterfile, and return the dictionary of the results
def read_pixel_masterfile(filename,usefloat=False):
	jd = 0
	array = []
	with open(filename) as f:
		for line in f:
			if jd == 0:
				jd = float(line.strip())
			else:
				val = line.strip().split()
				# Convert to integers
				if usefloat:
					val = [float(x) for x in val]
				else:
					val = [int(x) for x in val]
				array.append(val)
	return jd, array

# Read in all of the pixel masterfiles
def read_pixel_masterfiles(prefix,usefloat=False):
	jds = []
	arrays = []
	for i in range(1,500):
		print(prefix+format(i, '03d')+'.txt')
		try:
			jd, array = read_pixel_masterfile(prefix+format(i, '03d')+'.txt',usefloat=usefloat)
		except:
			break
		jds.append(jd)
		arrays.append(array)
	return jds, arrays

# Read in a pixel positions file, and return the dictionary of the results
def read_pixel_positions(filename):
	array = []
	with open(filename) as f:
		for line in f:
			if '*' not in line:
				val = line.strip().split()
				array.append([int(val[0]), float(val[1]), float(val[2])])
	return array

def read_tod_files(indir, prefix, numpixels, numfiles=50,quiet=True):
	# Read in the data
	jd = np.empty(0)
	az = np.empty(0)
	el = np.empty(0)
	data = np.empty((0,0,0))
	for i in range(0,numfiles):
		print(indir+'/'+prefix+'-'+format(i, '04d')+'.tod2')
		try:
			inputfits = fits.open(indir+'/'+prefix+'-'+format(i, '04d')+'.tod2')
		except:
			break
		# cols = inputfits[1].columns
		# col_names = cols.names
		ndata = len(inputfits[1].data.field(0)[0][:])
		nsamples = ndata//4
		if nsamples*4 != ndata:
			print('Oddity in ' + prefix+'-'+format(i, '04d')+'.tod2' + ' - ndata = ' + str(ndata) + ' is not dividable by 4, changing it to ' + str(nsamples*4))
		if i == 0:
			jd = inputfits[1].data.field(0)[0][:nsamples*4]
			az = inputfits[1].data.field(1)[0][:nsamples*4]
			el = inputfits[1].data.field(2)[0][:nsamples*4]
		else:
			jd = np.append(jd, inputfits[1].data.field(0)[0][:nsamples*4])
			az = np.append(az,inputfits[1].data.field(1)[0][:nsamples*4])
			el = np.append(el,inputfits[1].data.field(2)[0][:nsamples*4])
		rawdata = inputfits[1].data.field(3)[0][:nsamples*4*4*31]
		if np.shape(rawdata)[0] == ndata:
			rawdata = rawdata.reshape(ndata*124,order='C')
			rawdata = rawdata[:nsamples*4*4*31]
		data = np.append(data, rawdata)

		if not quiet:
			print(' Start time: ' + str(np.min(inputfits[1].data.field(0)[0][:nsamples*4])))
			print(' End time: ' + str(np.max(inputfits[1].data.field(0)[0][:nsamples*4])))
			print(' Duration: ' + str((np.max(inputfits[1].data.field(0)[0][:nsamples*4])-np.min(inputfits[1].data.field(0)[0][:nsamples*4]))*24*60*60) + ' seconds')
			print(' There are ' + str(ndata) + " datapoints")
			print(' Az range: ' + str(np.min(inputfits[1].data.field(1)[0][:nsamples*4])) + ' to ' + str(np.max(inputfits[1].data.field(1)[0][:nsamples*4])))
			print(' El range: ' + str(np.min(inputfits[1].data.field(2)[0][:nsamples*4])) + ' to ' + str(np.max(inputfits[1].data.field(2)[0][:nsamples*4])))
			print(' Raw data array is ' + str(np.shape(rawdata)))

	# Reshape the data
	ndata = len(jd)//4
	az = az.reshape(4,ndata,order='F')
	el = el.reshape(4,ndata,order='F')
	jd = jd.reshape(4,ndata,order='F')
	print(prefix)
	print(' Start time: ' + str(np.min(jd)))
	print(' End time: ' + str(np.max(jd)))
	print(' Duration: ' + str((np.max(jd)-np.min(jd))*24*60*60) + ' seconds')
	print(' There are ' + str(ndata) + " datapoints")
	print(' Az range: ' + str(np.min(az)) + ' to ' + str(np.max(az)))
	print(' El range: ' + str(np.min(el)) + ' to ' + str(np.max(el)))
	print(' JD shape: ' + str(np.shape(jd)))
	print(' Az shape: ' + str(np.shape(az)))
	print(' Data shape: ' + str(np.shape(data)))
	print(' New data length: ' + str(len(data)/(numpixels*4*4)))
	data = data.reshape(4, numpixels, 4, ndata, order='F')
	print(' New data shape: ' + str(np.shape(data)))
	return az, el, jd, data
