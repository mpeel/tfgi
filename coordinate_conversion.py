#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Sample script to convert Az/El to RA/Dec
# 
# Version history:
#
# 29-May-2019  M. Peel       Started

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, Angle
from astropy.time import Time
import astropy.units as u
import numpy as np
import time

# For astropy
telescope = EarthLocation(lat=28.300224*u.deg, lon=-16.510113*u.deg, height=2390*u.m)

# Set up the dataset
numsamples = 100000
times = 2458000.5 + np.arange(0,numsamples)*4.6e-8
az = 100.0 + np.arange(0,numsamples)*((200-100)/numsamples)
el = 30.0 + np.arange(0,numsamples)*((50-30)/numsamples)
print(times)
print(az)
print(el)

# This is the conversion function
def convert_azel_radec(az,el,timearr):
	return AltAz(az=az*u.deg,alt=el*u.deg,location=telescope,obstime=timearr).transform_to(ICRS)

# Start the timing
start = time.time()

# Do the conversion
timearr = Time(times, format='jd')
positions = convert_azel_radec(az,el,timearr)
print(positions)

# Print the runtime
end = time.time()
print(str(end-start) + ' seconds')

