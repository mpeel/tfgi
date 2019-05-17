#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Various functions related to tgi, separate from the class
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

def linfit(x, A, B):
	return A*x+B

def fit_kneefreq(freq, param):
	sigma, fknee, alpha = param
	return sigma**2 * (1 + (fknee / freq)**alpha)

def compute_residuals(param, data, freq):
	model = fit_kneefreq(freq, param)
	residual = np.log(data / model)
	return residual

def fit_skydip(el, param):
	return param[0] + param[1]/(np.sin(el*np.pi/180.0))# + param[2]*el

def compute_residuals_skydip(param, el, data):
	model = fit_skydip(el, param)
	residual = data - model
	return residual

def ensure_dir(f):
	os.makedirs(f, exist_ok=True)

def fit_gaussian(x, param):
    return param[0]*np.exp(-(x*x)/(2*param[1]*param[1]))+param[2]

def compute_residuals_gaussian(param, x, y):
	return y - fit_gaussian(x,param)
