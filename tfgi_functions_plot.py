#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Various functions related to plotting tgi data, separate from the class
# 
# Version history:
#
# 17-May-2019  M. Peel       Split from tgi.py

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
from tfgi_functions import *

def plot_tfgi_tod(data, outputname,formatstr='b.'):
	plt.plot(data,formatstr)
	plt.xlabel('Samples')
	plt.ylabel('Power')
	plt.savefig(outputname)
	plt.close()
	plt.clf()
	return

# Plot the TODs against a given set of vals, e.g. az or el.
def plot_tfgi_val_tod(val, data, outputname):
	plt.plot(val,data,'b.')
	plt.savefig(outputname)
	plt.close()
	plt.clf()
	return

# Plot a skydip with a fit
def plot_tfgi_skydip(el,data,outputname):
	# Fit the skydip
	params = [1,1,0]
	param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_skydip, params, args=(el, data),full_output=True)
	sigma_param_est = np.sqrt(np.diagonal(cov_x))
	# Now plot things
	plt.plot(el,data,'b.')
	mesg_fit = (
	r'$A={:5.3g}\pm{:3.2g}$'.format(
		param_est[0], sigma_param_est[0]) + ','
	r'$B={:5.3f}\pm{:3.2f}$'.format(
		param_est[1], sigma_param_est[1]) + ','
	r'     $C={:5.3f}\pm{:3.2f}$'.format(
		param_est[2], sigma_param_est[2]))
	plt.plot(el,fit_skydip(el,param_est),'g',label="Fit: " + mesg_fit)
	plt.legend(prop={'size':8})
	plt.savefig(outputname)
	plt.close()
	plt.clf()
	return param_est
