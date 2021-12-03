import astropy.units as u
from astropy.time import Time
from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz

# Take any ICRS sky coordinate
icrs = SkyCoord.from_name('crab')
print('RA = {pos.ra.deg:10.5f}, DEC = {pos.dec.deg:10.5f}'.format(pos=icrs))
# RA =   83.63308, DEC =   22.01450

# Convert to AltAz for some random observation time and location
# This assumes pressure is zero, i.e. no refraction
time = Time('2010-04-26', scale='tt')
location = EarthLocation(lat=28.300224*u.deg, lon=-16.510113*u.deg, height=2390*u.m)
altaz_frame = AltAz(obstime=time, location=location)
altaz = icrs.transform_to(altaz_frame)
print('AZ = {pos.az.deg:10.5f}, ALT = {pos.alt.deg:10.5f}'.format(pos=altaz))
# AZ =  351.88232, ALT =  -25.56281

# Convert back to ICRS to make sure round-tripping is OK
icrs2 = altaz.transform_to('icrs')
print('RA = {pos.ra.deg:10.5f}, DEC = {pos.dec.deg:10.5f}'.format(pos=icrs2))
# RA =   83.63308, DEC =   22.01450