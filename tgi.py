#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Various functions related to tgi
# 
# Version history:
#
# 15-Apr-2019  M. Peel       Started

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
from tgi_functions import *
from tgi_functions_read import *
import logging

def calc_beam_area(beam):
	return (np.pi * (beam*np.pi/180.0)**2)/(4.0*np.log(2.0))

def calc_Tsys(std,B,t):
	return std*np.sqrt(B*t)

class tgi:
	def __init__(self,indir="/net/nas/proyectos/quijote2",outdir='',pixelfileloc="/net/nas/proyectos/quijote2/etc/qt2_pixel_masterfile.",pixelposfileloc="/net/nas/proyectos/quijote2/etc/tgi_fgi_horn_positions_table.txt",polcalfileloc="/net/nas/proyectos/quijote2/etc/qt2_polcal."):
		self.numpixels = 31
		self.indir = indir
		self.outdir = outdir
		self.telescope = EarthLocation(lat=28.300224*u.deg, lon=-16.510113*u.deg, height=2390*u.m)
		self.jd_ref = 2456244.5 # Taken as the first day of the Commissioning (13/Nov/2012, 0.00h)
		self.nside = 512
		self.npix = hp.nside2npix(self.nside)
		self.pixeljds, self.pixelarrays = self.read_pixel_masterfiles(pixelfileloc)
		self.polcaljds, self.polcalarrays = self.read_pixel_masterfiles(polcalfileloc,usefloat=True)
		# print(self.polcaljds, self.polcalarrays)
		# exit()
		self.pixelpositions = self.read_pixel_positions(pixelposfileloc)
		self.apobserver = astroplan.Observer(latitude=28.300224*u.deg, longitude=-16.510113*u.deg, elevation=2390*u.m)
		# NB: Galactic doesn't work with hour angles
		self.coordsys = 1 # 0 = Galactic, 1 = ICRS.

		self.k = 1.380e-23 # Boltzman constant m2 kg s-2 K-1
		self.c = 2.9979e8

		# Roughly
		self.nu_tgi = 30e9
		self.B_tgi = 8e9
		self.nu_fgi = 40e9
		self.B_fgi = 8e9

	def calc_JytoK(self,beam,freq):
		return 1e-26*((self.c/freq)**2)/(2.0*self.k*calc_beam_area(beam))

	# Read in a single pixel masterfile, and return the dictionary of the results
	def read_pixel_masterfile(self, filename,usefloat=False):
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
	def read_pixel_masterfiles(self, prefix,usefloat=False):
		jds = []
		arrays = []
		for i in range(1,500):
			print(prefix+format(i, '03d')+'.txt')
			try:
				jd, array = self.read_pixel_masterfile(prefix+format(i, '03d')+'.txt',usefloat=usefloat)
			except:
				break
			jds.append(jd)
			arrays.append(array)
		# print(jds)
		# exit()
		return jds, arrays

	# Read in a pixel positions file, and return the dictionary of the results
	def read_pixel_positions(self, filename):
		array = []
		with open(filename) as f:
			for line in f:
				if '*' not in line:
					val = line.strip().split()
					array.append([int(val[0]), float(val[1]), float(val[2])])
		return array

	def get_pixel_info(self, jd, das):
		pixel_jd = [x for x in self.pixeljds if x <= jd][-1]
		# print(pixel_jd)
		pixel_id = self.pixeljds.index(pixel_jd)
		# print(pixel_id)
		for line in self.pixelarrays[pixel_id]:
			# print(line)
			if line[5] == das:
				if line[1] >= 40:
					tgi = 0
					fgi = 1
				else:
					tgi = 1
					fgi = 0
				pixelposition = self.pixelpositions[line[0]-1]
				# print(pixelposition)

				# See if we can get polcal info
				polcal = self.get_polcal_info(jd,das)

				return {'fp':line[0],'pixel':line[1],'fem':line[2],'bem':line[3],'das':line[5],'tgi':tgi,'fgi':fgi,'x_pos':pixelposition[1],'y_pos':pixelposition[2],'polcal':polcal}
				# return line
		return []

	def get_polcal_info(self, jd, das):
		try:
			pixel_jd = [x for x in self.polcaljds if x <= jd][-1]
		except:
			return []
		pixel_id = self.polcaljds.index(pixel_jd)
		returnvals = []
		for line in self.polcalarrays[pixel_id]:
			if line[0] == das:
				returnvals.append([line[1],line[2],line[3]])#,line[4],line[5],line[6],line[7]])
		return returnvals

	def read_tod(self, prefix, numfiles=50,quiet=True):
		return read_tod_files(self.indir, prefix, self.numpixels, numfiles=numfiles,quiet=quiet)

	# Note: this is currently slow and over-precise, see
	# https://github.com/astropy/astropy/pull/6068
	def calc_positions(self, az, el, jd):
		# Convert the az/el in the data to healpix pixel numbers
		time = Time(jd[0]+self.jd_ref, format='jd')
		position = AltAz(az=az[0]*u.deg,alt=el[0]*u.deg,location=self.telescope,obstime=time)
		if self.coordsys == 0:
			skypos = position.transform_to(Galactic)
			pa = []
		else:
			skypos = position.transform_to(ICRS)
			# pa = []
			pa = self.apobserver.parallactic_angle(time,skypos)
		return skypos,pa

	def calc_healpix_pixels(self, skypos):
		if self.coordsys == 0:
			healpix_pixel = hp.ang2pix(self.nside, (np.pi/2)-Angle(skypos.b).radian, Angle(skypos.l).radian)
			pos = (np.median(Angle(skypos.l).degree),np.median((Angle(skypos.b).degree)))
		else:
			print(skypos)
			healpix_pixel = hp.ang2pix(self.nside, (np.pi/2)-Angle(skypos.dec).radian, Angle(skypos.ra).radian)
			pos = (np.median(Angle(skypos.ra).degree),np.median((Angle(skypos.dec).degree)))
		return healpix_pixel, pos

	def plot_tod(self, data, outputname,formatstr='b.'):
		plt.plot(data,formatstr)
		plt.xlabel('Samples')
		plt.ylabel('Power')
		plt.savefig(outputname)
		plt.close()
		plt.clf()
		return

	# Plot the TODs against a given set of vals, e.g. az or el.
	def plot_val_tod(self, val, data, outputname):
		plt.plot(val,data,'b.')
		plt.savefig(outputname)
		plt.close()
		plt.clf()
		return

	# Plot a skydip with a fit
	def plot_skydip(self,el,data,outputname):
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

	def plot_fft(self, data, outputname,samplerate=1000,numsmoothbins=50):
		fft_w = np.fft.rfft(data)
		fft_f = np.fft.rfftfreq(len(data), 1/samplerate)
		fft_spectrum = np.abs(fft_w)
		#fft_spectrum = np.abs(fft_w)**2
		# plt.plot(fft_f, fft_spectrum)
		# plt.xscale("log")
		# plt.yscale("log")
		# plt.savefig(self.outdir+'/'+prefix+'/plot_fft_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
		# plt.clf()
		# plt.plot(fft_f, fft_spectrum)
		# plt.yscale("log")
		# plt.savefig(self.outdir+'/'+prefix+'/plot_fft_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_lin.png')
		# plt.clf()

		# Do a smoothed version
		bins = np.linspace(np.log(min(fft_f[1:-1])),np.log(max(fft_f[1:-1])),numsmoothbins,endpoint=False)
		values = np.zeros(numsmoothbins)
		values2 = np.zeros(numsmoothbins)
		matches = np.digitize(np.log(fft_f),bins)
		binmask = np.ones(numsmoothbins)
		for i in range(0,numsmoothbins):
			values[i] = np.nanmean(np.log(fft_spectrum[matches == i]))
			values2[i] = np.nanmean(fft_spectrum[matches == i])
			if np.isnan(values[i]) or np.isnan(values2[i]):
				binmask[i] = 0
		binmask[0] = 0
		binmask[-1] = 0
		binmask[-2] = 0
		# Fit to get the knee frequency
		# Using modified code from http://pchanial.github.io/python-for-data-scientists/body.html (by Jm. Colley)
		# params = np.array([np.sqrt(np.median(fft_spectrum))/5, 0.01, 1.0])
		params = np.array([(values[binmask==1])[-2], 0.01, 1.0])
		#print(params)
		# print np.exp(values[0,5:])
		# print np.exp(bins[5:])
		#print(values[binmask==1])
		#print(values2[binmask==1])
		# param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals, params, args=(values2[binmask==1], bins[binmask==1]),full_output=True)
		param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals, params, args=(np.exp(values[binmask==1]), np.exp(bins[binmask==1])),full_output=True)
		#print(param_est)
		#print(cov_x)
		#print(mesg_result)
		#print(ret_value)
		# exit()

		# Plot the spectrum
		plt.xlabel('Frequency (Hz)')
		plt.ylabel('Power')
		plt.xscale('log')
		plt.yscale('log')
		ymax = np.max(fft_spectrum[1:])*1.5
		ymin = np.min(fft_spectrum[1:])
		plt.ylim(ymin=ymin,ymax=ymax)
		plt.plot(fft_f, fft_spectrum,label='Data')
		plt.plot(np.exp(bins[1:]), np.exp(values[1:]), 'r', label='Mean (in log)')
		plt.plot(np.exp(bins[1:]), values2[1:], 'm', label='Mean (in linear)')
		try:
			sigma_param_est = np.sqrt(np.diagonal(cov_x))
			mesg_fit = (
			r'$\sigma={:5.3g}\pm{:3.2g}$'.format(
				param_est[0], sigma_param_est[0]) + ','
			r'$f_{{\rm knee}}={:5.3f}\pm{:3.2f}$'.format(
				param_est[1], sigma_param_est[1]) + ','
			r'     $\alpha={:5.3f}\pm{:3.2f}$'.format(
				param_est[2], sigma_param_est[2]))
			plt.plot(fft_f, fit_kneefreq(fft_f, param_est),label="Fit: " + mesg_fit)
		except:
			mesg_fit = (
			r'$\sigma={:5.3g}$'.format(
				param_est[0]) + ','
			r'$f_{{\rm knee}}={:5.3f}$'.format(
				param_est[1]) + ','
			r'     $\alpha={:5.3f}$'.format(
				param_est[2]))
			plt.plot(fft_f, fit_kneefreq(fft_f, param_est),label="Fit" + mesg_fit)
		plt.legend(prop={'size':8})
		plt.savefig(outputname)
		plt.close()
		plt.clf()

		return param_est, sigma_param_est

	def subtractbaseline(self, data, option=0, navg=1500):
		# Crude baseline removal
		npoints = len(data)
		# print(npoints)
		for i in range(0,npoints//navg):
			start = i*navg
			end = ((i+1)*navg)
			# print('start:' + str(start) + ', end: ' + str(end))
			if option == 0:
				data[start:end] = data[start:end] - np.median(data[start:end])
			else:
				xline = range(0,len(data[start:end]))
				if len(xline) != 0:
					A,B=curve_fit(linfit,xline,data[start:end])[0]
					# print(A,B)
					data[start:end] = data[start:end] - (A*xline+B)

		if option == 0:
			data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - np.median(data[(npoints//navg)*navg-1:])
		else:
			xline = range(0,len(data[(npoints//navg)*navg-1:]))
			if len(xline) != 0:
				A,B=curve_fit(linfit,xline,data[(npoints//navg)*navg-1:])[0]
				# print(A,B)
				data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - (A*xline+B)

		return data

	def converttopol(self, data, ordering=[0,1,2,3],pa=[],polangle=0,polcorr=[1.0,1.0,1.0,1.0]):
		newdata = np.zeros(np.shape(data))

		# This didn't work, so isn't currently used.
		# data[ordering[0]] = data[ordering[0]] * polcorr[0]
		# data[ordering[1]] = data[ordering[1]] * polcorr[1]
		# data[ordering[2]] = data[ordering[2]] * polcorr[2]
		# data[ordering[3]] = data[ordering[3]] * polcorr[3]

		newdata[0] = (data[ordering[0],:] + data[ordering[1],:])# / 2.0
		newdata[1] = (data[ordering[0],:] - data[ordering[1],:])# / 2.0
		newdata[2] = (data[ordering[2],:] - data[ordering[3],:])# / 2.0
		newdata[3] = (data[ordering[2],:] + data[ordering[3],:])# / 2.0

		if len(pa) > 0:
			ang = 2.0*pa # Already in radians
			Q = newdata[1,:]*np.cos(ang) + newdata[2,:]*np.sin(ang)
			U = -newdata[1,:]*np.sin(ang) + newdata[2,:]*np.cos(ang) 
			newdata[1] = Q.copy()
			newdata[2] = U.copy()

		# print(polangle)
		if polangle != 0:
			Q = newdata[1,:]*np.cos(2.0*polangle*np.pi/180.0) + newdata[2,:]*np.sin(2.0*polangle*np.pi/180.0)
			U = -newdata[1,:]*np.sin(2.0*polangle*np.pi/180.0) + newdata[2,:]*np.cos(2.0*polangle*np.pi/180.0)
			newdata[1] = Q.copy()
			newdata[2] = U.copy()

		return newdata


	def analyse_tod(self, prefix,pixelrange=range(0,31),detrange=range(0,4),phaserange=range(0,4),plotlimit=0.0,quiet=False,dofft=False,plottods=True,plotmap=True,dopol=False,plotcombination=False,numfiles=50):
		print(self.outdir+'/'+prefix)
		ensure_dir(self.outdir+'/'+prefix)

		if 'DIP' in prefix:
			# We have a sky dip. Run the routine to analyse that rather than this routine.
			self.analyse_skydip(prefix,pixelrange=pixelrange,detrange=detrange,phaserange=phaserange,quiet=quiet,plottods=plottods,dopol=dopol)
			return

		# Read in the data
		az, el, jd, data = self.read_tod(prefix,numfiles=numfiles,quiet=quiet)

		# Calculate the Galactic Healpix pixels and the central position from az/el
		skypos,pa = self.calc_positions(az, el, jd)
		healpix_pixel, centralpos = self.calc_healpix_pixels(skypos)

		self.plot_tod(skypos.ra,self.outdir+'/'+prefix+'/plot_ra.png')
		self.plot_tod(skypos.dec,self.outdir+'/'+prefix+'/plot_dec.png')
		self.plot_tod(pa,self.outdir+'/'+prefix+'/plot_pa.png')
		# exit()

		# Make maps for each pixel, detector, phase
		for pix in pixelrange:
			# Get the pixel info
			pixinfo = self.get_pixel_info(jd[0][0]+self.jd_ref,pix+1)
			print(pixinfo)
			if pixinfo == []:
				# We don't have a good pixel, skip it
				continue
			if pixinfo['pixel'] <= 0:
				# Pixel isn't a pixel
				continue

			for det in detrange:

				for j in phaserange:
					print(j)
					if plottods:
						# Plot some tods
						self.plot_tod(data[det][pix][j][10:-1500], self.outdir+'/'+prefix+'/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_pre.png')
						self.plot_tod(data[det][pix][j][10:5000],self.outdir+'/'+prefix+'/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_pre_zoom.png')
					# Do an FFT. Warning, this can slow things down quite a bit.
					if dofft:
						param_est, sigma_param_est = self.plot_fft(data[det][pix][j][0:10000],self.outdir+'/'+prefix+'/plot_fft_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_fit.png',samplerate=1000)
						print(str(param_est[0]) + " " + str(sigma_param_est[0]) + " " + str(param_est[1]) + " " + str(sigma_param_est[1]) + " " + str(param_est[2]) + " " + str(sigma_param_est[2]))
					data[det][pix][j] = self.subtractbaseline(data[det][pix][j])
					if plottods:
						self.plot_tod(data[det][pix][j][10:-1500], self.outdir+'/'+prefix+'/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						self.plot_tod(data[det][pix][j][10:5000], self.outdir+'/'+prefix+'/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_zoom.png')

				# Do we want to change from phase to I1,Q,U,I2?
				if dopol:
					if pixinfo['tgi'] == 1:
						# if det == 0:
						ordering = [0,1,2,3]
						# elif det == 1:
						# 	ordering = [1,2,3,0]
						# elif det == 2:
						# 	ordering = [2,3,0,1]
						# elif det == 3:
						# 	ordering = [3,0,1,2]
					else:
						# if det == 0:
						ordering = [0,2,1,3]
						# elif det == 1:
						# 	ordering = [2,1,3,0]
						# elif det == 2:
						# 	ordering = [1,3,0,2]
						# elif det == 3:
						# 	ordering = [3,0,2,1]
					# If we have polcal data, use it
					print(pixinfo['polcal'])
					print(pixinfo['polcal'][det])
					if pixinfo['polcal'] != []:
						# polcorr = [pixinfo['polcal'][det][3],pixinfo['polcal'][det][4],pixinfo['polcal'][det][5],pixinfo['polcal'][det][6]]
						data[det][pix] = self.converttopol(data[det][pix],ordering=ordering,pa=pa,polangle=pixinfo['polcal'][det][1])#,polcorr=polcorr)
					else:
						data[det][pix] = self.converttopol(data[det][pix],ordering=ordering,pa=pa)

				for j in phaserange:
					print('Pixel ' + str(pix+1) + ', detector ' + str(det+1) + ', phase ' + str(j+1))

					if plottods:
						self.plot_tod(data[det][pix][j][10:-1500], self.outdir+'/'+prefix+'/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_cal.png')
						self.plot_tod(data[det][pix][j][10:5000], self.outdir+'/'+prefix+'/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_cal_zoom.png')

					skymap = np.zeros(self.npix, dtype=np.float)
					hitmap = np.zeros(self.npix, dtype=np.float)
					for i in range(0,len(healpix_pixel)):
						skymap[healpix_pixel[i]] = skymap[healpix_pixel[i]] + data[det][pix][j][i]
						hitmap[healpix_pixel[i]] = hitmap[healpix_pixel[i]] + 1
					for i in range(0,len(skymap)):
						if hitmap[i] >= 1:
							skymap[i] = skymap[i]/hitmap[i]
						else:
							skymap[i] = hp.pixelfunc.UNSEEN

					# Get the maximum value in the map
					maxval = np.max(skymap)
					# Get a sample of data from the start of the measurement
					std = np.std(data[det][pix][j][0:4000])
					print(maxval)
					print(std)
					if pixinfo['tgi']:
						conv = self.calc_JytoK(0.2,self.nu_tgi)
						flux = 344.0
					else:
						conv = self.calc_JytoK(0.2,self.nu_fgi)
						flux = 318.0

					estimate = std*(flux/maxval)/np.sqrt(1000.0)
					print('In Jy/sec:' + str(estimate))
					print('In K/sec:' + str(estimate*conv))
					print('System temperature:' + str(calc_Tsys(estimate*conv,8e9,1/1000)))

					hp.write_map(self.outdir+'/'+prefix+'/skymap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.fits',skymap,overwrite=True)
					hp.write_map(self.outdir+'/'+prefix+'/hitmap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.fits',hitmap,overwrite=True)

					if plotmap:
						hp.mollview(skymap)
						plt.savefig(self.outdir+'/'+prefix+'/skymap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						plt.close()
						plt.clf()
					if plotmap:
						hp.mollview(hitmap)
						plt.savefig(self.outdir+'/'+prefix+'/hitmap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						hp.gnomview(skymap,rot=centralpos,reso=5)
						plt.savefig(self.outdir+'/'+prefix+'/skymap_zoom_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						plt.close()
						plt.clf()
					if plotmap and plotlimit != 0.0:
						hp.mollview(skymap,min=-plotlimit,max=plotlimit)
						plt.savefig(self.outdir+'/'+prefix+'/skymap_cut_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						plt.close()
						plt.clf()
						hp.gnomview(skymap,rot=centralpos,reso=5,min=-plotlimit,max=plotlimit)
						plt.savefig(self.outdir+'/'+prefix+'/skymap_zoom_cut_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						plt.close()
						plt.clf()


		# Do some combined plots if requested
		if plotcombination:
			for i in range(1,5):
				if i == 2 or i == 3:
					plotlimit2 = plotlimit
				else:
					plotlimit2 = 0.0
				for k in range(1,5):
					filelist = []
					hitlist = []
					for pix in range(1,self.numpixels):
						filelist.append(self.outdir+'/'+prefix+'/skymap_'+str(pix)+'_'+str(k)+'_'+str(i)+'.fits')
						hitlist.append(self.outdir+'/'+prefix+'/hitmap_'+str(pix)+'_'+str(k)+'_'+str(i)+'.fits')
					print(filelist)
					run.combine_sky_maps(filelist,hitlist,self.outdir+'/'+prefix+'/combined_'+str(i)+'_'+str(k),centralpos=centralpos,plotlimit=plotlimit2)

			for i in range(1,5):
				run.calc_P_angle_skymaps(self.outdir+'/'+prefix+'/combined_1_'+str(i)+'_skymap.fits',self.outdir+'/'+prefix+'/combined_2_'+str(i)+'_skymap.fits',self.outdir+'/'+prefix+'/combined_3_'+str(i)+'_skymap.fits',self.outdir+'/'+prefix+'/combined_pol_'+str(i),centralpos=centralpos)


		# # Plot a boxcar average
		# numboxcar = 101
		# boxcar = signal.boxcar(numboxcar)
		# boxcar = boxcar/np.sum(boxcar)

		# print(np.shape(data))
		# for i in range(0,31):
		# 	print(data[i])
		# 	print (len(data[i]))
		# 	fig = plt.figure(figsize=(20,10))
		# 	for j in range(0,4):
		# 		toplot = signal.convolve(data[i][j],boxcar)
		# 		print(len(toplot))
		# 		plt.plot(toplot[numboxcar:ndata-numboxcar]-np.median(toplot[numboxcar:ndata-numboxcar]))
		# 	plt.savefig('plot_boxcar_'+str(i+1)+'_tod.png')
		# 	plt.clf()

		return

	def stack_maps_tod(self, prefix,schedules,pixelrange=range(0,31),detrange=range(0,4),phaserange=range(0,4),plotlimit=0.0,quiet=False,dofft=False,dopol=False,numfiles=50):
		print(self.outdir+'/'+prefix)
		ensure_dir(self.outdir+'/'+prefix)

		skymap = np.zeros((len(pixelrange),self.npix), dtype=np.float)
		hitmap = np.zeros((len(pixelrange),self.npix), dtype=np.float)
		centralpos = (0,0)
		for schedule in schedules:
			# Read in the data
			az, el, jd, data = self.read_tod(schedule,numfiles=numfiles,quiet=quiet)

			# Calculate the Galactic Healpix pixels and the central position from az/el
			healpix_pixel, centralpos = self.calc_positions(az, el, jd)

			# Make maps for each pixel, detector, phase
			pixnum=-1
			for pix in pixelrange:
				pixnum = pixnum+1
				for det in detrange:
					# Do we want to change from phase to I1,Q,U,I2?
					if dopol:
						# Different for TGI and FGI...
						if pix > 10:
							option = 1
						else:
							option = 0
						# data[det][pix] = self.converttopol(data[det][pix],option=option)

					for j in phaserange:
						data[det][pix][j] = self.subtractbaseline(data[det][pix][j])
						for i in range(0,len(healpix_pixel)):
							skymap[pixnum][healpix_pixel[i]] = skymap[pixnum][healpix_pixel[i]] + data[det][pix][j][i]
							hitmap[pixnum][healpix_pixel[i]] = hitmap[pixnum][healpix_pixel[i]] + 1

		# We're done making the maps, now normalise them and write it out
		pixnum=-1
		for pix in pixelrange:
			pixnum = pixnum+1
			for i in range(0,len(skymap[pixnum])):
				if hitmap[pixnum][i] >= 1:
					skymap[pixnum][i] = skymap[pixnum][i]/hitmap[pixnum][i]
				else:
					skymap[pixnum][i] = hp.pixelfunc.UNSEEN

			hp.write_map(self.outdir+'/'+prefix+'/skymap_'+str(pix+1)+'.fits',skymap[pixnum],overwrite=True)
			hp.mollview(skymap[pixnum])
			plt.savefig(self.outdir+'/'+prefix+'/skymap_'+str(pix+1)+'.png')
			plt.close()
			plt.clf()
			hp.write_map(self.outdir+'/'+prefix+'/hitmap_'+str(pix+1)+'.fits',hitmap[pixnum],overwrite=True)
			hp.mollview(hitmap[pixnum])
			plt.savefig(self.outdir+'/'+prefix+'/hitmap_'+str(pix+1)+'.png')
			hp.gnomview(skymap[pixnum],rot=centralpos,reso=5)
			plt.savefig(self.outdir+'/'+prefix+'/skymap_zoom_'+str(pix+1)+'.png')
			plt.close()
			plt.clf()
			if plotlimit != 0.0:
				hp.mollview(skymap[pixnum],min=-plotlimit,max=plotlimit)
				plt.savefig(self.outdir+'/'+prefix+'/skymap_cut_'+str(pix+1)+'.png')
				plt.close()
				plt.clf()
				hp.gnomview(skymap[pixnum],rot=centralpos,reso=5,min=-plotlimit,max=plotlimit)
				plt.savefig(self.outdir+'/'+prefix+'/skymap_zoom_cut_'+str(pix+1)+'.png')
				plt.close()
				plt.clf()
		return

	def combine_sky_maps(self,skymaps,hitmaps,outputname,centralpos=(0,0),plotlimit=0.0):

		skymap = np.zeros(self.npix, dtype=np.float)
		hitmap = np.zeros(self.npix, dtype=np.float)

		nummaps = len(skymaps)
		for i in range(0,nummaps):
			try:
				inputmap = hp.read_map(skymaps[i])
				inputhitmap = hp.read_map(hitmaps[i])
			except:
				continue
			for j in range(0,self.npix):
				if inputhitmap[j] > 0:
					skymap[j] = skymap[j] + inputmap[j]*inputhitmap[j]
					hitmap[j] = hitmap[j] + inputhitmap[j]

		# We now have a combined map, time to normalise it
		for i in range(0,self.npix):
			if hitmap[i] >= 1:
				skymap[i] = skymap[i]/hitmap[i]
			else:
				skymap[i] = hp.pixelfunc.UNSEEN

		hp.write_map(outputname+'_skymap.fits',skymap,overwrite=True)
		hp.mollview(skymap)
		plt.savefig(outputname+'_skymap.png')
		plt.close()
		plt.clf()
		hp.write_map(outputname+'_hitmap.fits',hitmap,overwrite=True)
		hp.mollview(hitmap)
		plt.savefig(outputname+'_hitmap.png')
		if plotlimit != 0.0:
			hp.gnomview(skymap,rot=centralpos,reso=5,min=-plotlimit,max=plotlimit)
		else:
			hp.gnomview(skymap,rot=centralpos,reso=5)
		plt.savefig(outputname+'_zoom.png')
		plt.close()
		plt.clf()

		return

	def calc_P_angle_skymaps(self,I_filename,Q_filename,U_filename,outputname,centralpos=(0,0),plotlimit=0.0):
		print(I_filename)
		I = hp.read_map(I_filename)
		Q = hp.read_map(Q_filename)
		U = hp.read_map(U_filename)
		P = np.sqrt(Q**2+U**2)
		ang = (0.5*np.arctan2(U,Q)* 180 / np.pi)
		frac = P/I

		P[I == hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
		ang[I == hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
		ang[P < 0.1*np.max(P)] = hp.pixelfunc.UNSEEN
		frac[I == hp.pixelfunc.UNSEEN] = hp.pixelfunc.UNSEEN
		frac[P < 0.1*np.max(P)] = hp.pixelfunc.UNSEEN

		hp.write_map(outputname+'_P.fits',P,overwrite=True)
		hp.mollview(P)
		plt.savefig(outputname+'_P.png')
		plt.close()
		plt.clf()
		hp.write_map(outputname+'_ang.fits',ang,overwrite=True)
		hp.mollview(ang)
		plt.savefig(outputname+'_ang.png')
		plt.close()
		plt.clf()
		hp.write_map(outputname+'_Pfrac.fits',frac,overwrite=True)
		hp.mollview(frac,min=0,max=10)
		plt.savefig(outputname+'_Pfrac.png')
		plt.close()
		plt.clf()
		if plotlimit != 0.0:
			hp.gnomview(P,rot=centralpos,reso=5,min=-plotlimit,max=plotlimit)
		else:
			hp.gnomview(P,rot=centralpos,reso=5)
		plt.savefig(outputname+'_P_zoom.png')
		plt.close()
		plt.clf()
		hp.gnomview(ang,rot=centralpos,reso=5,cmap=plt.get_cmap('hsv'))
		plt.savefig(outputname+'_ang_zoom.png')
		plt.close()
		plt.clf()
		hp.gnomview(frac,rot=centralpos,reso=5,min=0)#,max=0.1)
		plt.savefig(outputname+'_Pfrac_zoom.png')
		plt.close()
		plt.clf()

		return		

	def examine_source(self,skymaps,hitmaps,outputname,sourcepos=(0,0),plotlimit=0.0):

		skymap = np.zeros(self.npix, dtype=np.float)
		hitmap = np.zeros(self.npix, dtype=np.float)

		nummaps = len(skymaps)
		for i in range(0,nummaps):
			try:
				inputmap = hp.read_map(skymaps[i])
				inputhitmap = hp.read_map(hitmaps[i])
			except:
				continue

			if sourcepos == (0,0):
				# We have a map but no position, let's find the peak position
				skymax = np.max(inputmap)
				pixelpos = np.where(inputmap==skymax)[0][0]
				print(pixelpos)
				sourcepos = hp.pixelfunc.pix2ang(self.nside,pixelpos,lonlat=True)

			print(sourcepos)
			sourcecoord = SkyCoord(sourcepos[0]*u.deg, sourcepos[1]*u.deg, frame=Galactic)
			x_val = []
			y_val = []
			for i in range(0,self.npix):
				if inputmap[i] != hp.UNSEEN:
					if inputmap[i] != 0.0:
						pixelpos = hp.pixelfunc.pix2ang(self.nside,i,lonlat=True)
						pixelcoord = SkyCoord(pixelpos[0]*u.deg, pixelpos[1]*u.deg, frame=Galactic)
						sep=pixelcoord.separation(sourcecoord)
						# print(sep.degree,inputmap[i])
						if sep.degree < 2.0:
							x_val.append(sep.degree)
							y_val.append(inputmap[i])

			x_val = np.array(x_val)
			y_val = np.array(y_val)

			plt.plot(x_val,y_val,'b.')

			params = [10.0,1.0,0.0]
			param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_gaussian, params, args=(x_val, y_val),full_output=True)
			sigma_param_est = np.sqrt(np.diagonal(cov_x))
			mesg_fit = (
			r'$A={:5.3f}\pm{:3.2f}$'.format(
				param_est[0], sigma_param_est[0]) + ','
			r'$B={:5.3f}\pm{:3.2f}$'.format(
				param_est[1], sigma_param_est[1]) + ','
			r'$C={:5.3f}\pm{:3.2f}$'.format(
				param_est[2], sigma_param_est[2]))
			toplot = np.arange(0,2.0,0.01)
			plt.plot(toplot,fit_gaussian(toplot,param_est),label=mesg_fit)
			plt.xlabel('Distance')
			plt.ylabel('Power')
			plt.yscale('log')
			plt.ylim((1e-4,skymax))

			plt.legend(prop={'size':8})
			plt.savefig('test.pdf')
			plt.close()
			plt.clf()

		return

	def analyse_skydip(self, prefix,pixelrange=range(0,31),detrange=range(0,4),phaserange=range(0,4),plotlimit=0.0,quiet=False,dofft=False,plottods=True,dopol=False,numfiles=50,minel=35.0,maxel=85.0,numelbins=100,plotindividual=False):
		print(self.outdir+'/'+prefix)
		ensure_dir(self.outdir+'/'+prefix)

		atm_tgi = 5.0
		atm_fgi = 10.0

		# Start a log file
		logging.basicConfig(filename=self.outdir+'/'+prefix+'/log_skydip.log',level=logging.DEBUG)
		logging.info('Assuming atmosphere of ' + str(atm_tgi) + 'K for TGI and ' + str(atm_fgi) + 'K for FGI')
		# Read in the data
		az, el, jd, data = self.read_tod(prefix,numfiles=numfiles,quiet=quiet)

		# Create a mask of the different elevation dips
		tempmask = el[0].copy()
		tempmask[tempmask < minel] = 0.
		tempmask[tempmask > maxel] = 0.
		self.plot_tod(tempmask,self.outdir+'/'+prefix+'/mask.pdf')
		skydip_mask = tempmask.copy()
		# Find out whether we're going up or down in elevation
		for i in range(1,len(skydip_mask)):
			if tempmask[i] != 0.:
				if tempmask[i-1] < tempmask[i]:
					skydip_mask[i] = -1
				else:
					skydip_mask[i] = 1
		# To catch noisy bits in transitions, require that sets of testlen all have to have the same value
		testlen = 10
		for i in range(0,len((skydip_mask/testlen)-testlen)):
			if np.sum(np.abs(skydip_mask[i*testlen:i*testlen+testlen])) != testlen:
				skydip_mask[i*testlen:i*testlen+testlen] = 0
		# Count how many different sections we have, and update the mask so we can extract them.
		num_skydips = 1
		notzero = 0
		for i in range(0,len(skydip_mask)):
			if skydip_mask[i] != 0:
				skydip_mask[i] = skydip_mask[i] * num_skydips
				notzero = 1
			else:
				if notzero == 1:
					num_skydips = num_skydips+1
					notzero = 0
		# print(num_skydips)
		self.plot_tod(skydip_mask,self.outdir+'/'+prefix+'/mask_bit.pdf',formatstr='b')

		# Make maps for each pixel, detector, phase
		for pix in pixelrange:
			# Get the pixel info
			pixinfo = self.get_pixel_info(jd[0][0]+self.jd_ref,pix+1)
			print(pixinfo)
			if pixinfo == []:
				# We don't have a good pixel, skip it
				continue
			if pixinfo['pixel'] <= 0:
				# Pixel isn't a pixel
				continue

			for det in detrange:
				# Do we want to change from phase to I1,Q,U,I2?
				if dopol:
					if pixinfo['tgi'] == 1:
						ordering = [0,1,2,3]
					else:
						ordering = [0,2,1,3]
					data[det][pix] = self.converttopol(data[det][pix],ordering=ordering)

				for j in phaserange:

					if plottods:
						# Do a plot of az vs. tod
						self.plot_tod(data[det][pix][j],self.outdir+'/'+prefix+'/tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1))
						self.plot_val_tod(az[det],data[det][pix][j],self.outdir+'/'+prefix+'/az_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1))
						self.plot_val_tod(el[det],data[det][pix][j],self.outdir+'/'+prefix+'/el_tod_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1))
						self.plot_val_tod(az[det],el[det],self.outdir+'/'+prefix+'/az_el_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1))

					# Plot the individual skydips
					elbins = np.zeros((5,num_skydips,numelbins))
					for i in range(1,num_skydips):

						elstep = (maxel-minel)/float(numelbins)
						stepmask = np.zeros(len(skydip_mask))
						stepmask[skydip_mask == i] = 1
						stepmask[skydip_mask == -i] = 1
						for k in range(0,numelbins):
							elstepmin = minel+k*elstep
							elstepmax = minel+(k+1)*elstep
							elbins[0][i][k] = (elstepmin+elstepmax)/2.0
							elmask = np.zeros(len(skydip_mask))
							elmask[el[det] > elstepmin] = elmask[el[det] > elstepmin] + 0.5
							elmask[el[det] < elstepmax] = elmask[el[det] < elstepmax] + 0.5
							elmask[elmask < 1] = 0
							elbins[1][i][k] = np.mean(data[det][pix][j][stepmask*elmask==1])
							# elbins[2][i][k] = np.sum(data[det][pix][j][stepmask*elmask==1])

						# Before doing the std, subtract out a fit
						params = [1,1]#,0]
						param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_skydip, params, args=(elbins[0][i], elbins[1][i]),full_output=True)

						for k in range(0,numelbins):
							elstepmin = minel+k*elstep
							elstepmax = minel+(k+1)*elstep
							elmask = np.zeros(len(skydip_mask))
							elmask[el[det] > elstepmin] = elmask[el[det] > elstepmin] + 0.5
							elmask[el[det] < elstepmax] = elmask[el[det] < elstepmax] + 0.5
							elmask[elmask < 1] = 0

							tempdata = data[det][pix][j][stepmask*elmask==1] - fit_skydip(el[det][stepmask*elmask==1],param_est)
							elbins[2][i][k] = np.mean(tempdata)
							elbins[3][i][k] = np.std(tempdata)
							elbins[4][i][k] = fit_skydip(elbins[0][i][k],param_est)


						if plotindividual:
							self.plot_skydip(elbins[0][i],elbins[1][i],self.outdir+'/'+prefix+'/average_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_'+str(i))
							# self.plot_skydip(elbins[0][i],elbins[2][i],self.outdir+'/'+prefix+'/average_sub_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_'+str(i))
							self.plot_skydip(elbins[0][i],elbins[3][i],self.outdir+'/'+prefix+'/std_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_'+str(i))
							# self.plot_skydip(elbins[0][i],elbins[4][i],self.outdir+'/'+prefix+'/std_average_fit_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_'+str(i))

					# Combine the sky dips into one plot
					for i in range(1,num_skydips):
						if i in skydip_mask:
							plt.plot(elbins[0][i],elbins[1][i],'b')
						else:
							plt.plot(elbins[0][i],elbins[1][i],'r')
					plt.xlabel('Elevation')
					plt.ylabel('Power')
					plt.savefig(self.outdir+'/'+prefix+'/average_combine_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.pdf')
					plt.close()
					plt.clf()

					# Do the same for the standard deviations, with some fitting
					tofit_x_std = []
					tofit_y_std = []
					for i in range(1,num_skydips):
						if i in skydip_mask:
							plt.plot(elbins[0][i],elbins[3][i],'b')
						else:
							plt.plot(elbins[0][i],elbins[3][i],'r')
						tofit_x_std.append(elbins[0][i])
						tofit_y_std.append(elbins[3][i])
					tofit_x_std = np.array(tofit_x_std)
					tofit_y_std = np.array(tofit_y_std)
					params = [1,1]#,0]
					param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_skydip, params, args=(tofit_x_std.flatten(), tofit_y_std.flatten()),full_output=True)
					sigma_param_est = np.sqrt(np.diagonal(cov_x))
					skydip_fitted = fit_skydip(elbins[0][1],param_est)
					mesg_fit = (
					r'$A={:5.3e}\pm{:3.2e}$'.format(
						param_est[0], sigma_param_est[0]) + ','
					r'$B={:5.3e}\pm{:3.2e}$'.format(
						param_est[1], sigma_param_est[1]))# + ','
					# r'     $C={:5.3e}\pm{:3.2e}$'.format(
						# param_est[2], sigma_param_est[2]))
					plt.plot(elbins[0][1],skydip_fitted,'g',label="Fit: " + mesg_fit)

					if pixinfo['tgi'] == 1:
						systemp = (param_est[0]/param_est[1])*atm_tgi
					else:
						systemp = (param_est[0]/param_est[1])*atm_fgi

					print('DAS '+str(pix) + ' Tsys estimate: ' + str(systemp))

					plt.xlabel('Elevation')
					plt.ylabel('Power')
					plt.legend(prop={'size':8})
					plt.savefig(self.outdir+'/'+prefix+'/std_combine_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.pdf')
					plt.close()
					plt.clf()

					# Do a fit to the data, and plot that along with the normalized data
					tofit_x = []
					tofit_y = []
					for i in range(1,num_skydips):
						newmax=1.0
						newmin=0.0
						max=np.max(elbins[1][i])
						min=np.min(elbins[1][i])
						toplot = (newmax-newmin)/(max-min)*(elbins[1][i]-max)+newmax
						tofit_x.append(elbins[0][i][:])
						tofit_y.append(toplot[:])
						if i in skydip_mask:
							plt.plot(elbins[0][i],toplot,'b',alpha=0.5)
						else:
							plt.plot(elbins[0][i],toplot,'r',alpha=0.5)
					tofit_x = np.array(tofit_x)
					tofit_y = np.array(tofit_y)
					params = [1,1]#,0]
					param_est, cov_x, infodict, mesg_result, ret_value = optimize.leastsq(compute_residuals_skydip, params, args=(tofit_x.flatten(), tofit_y.flatten()),full_output=True)
					sigma_param_est = np.sqrt(np.diagonal(cov_x))
					skydip_fitted = fit_skydip(elbins[0][1],param_est)
					mesg_fit = (
					r'$A={:5.3e}\pm{:3.2e}$'.format(
						param_est[0], sigma_param_est[0]) + ','
					r'$B={:5.3e}\pm{:3.2e}$'.format(
						param_est[1], sigma_param_est[1]))# + ','
					# r'     $C={:5.3e}\pm{:3.2e}$'.format(
						# param_est[2], sigma_param_est[2]))
					plt.plot(elbins[0][1],skydip_fitted,'g',label="Fit: " + mesg_fit)
					plt.xlabel('Elevation')
					plt.ylabel('Power')
					plt.legend(prop={'size':8})
					plt.savefig(self.outdir+'/'+prefix+'/average_combine_norm_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.pdf')
					plt.close()
					plt.clf()


					# Plot the residuals after subtracting the best fit
					for i in range(1,num_skydips):
						newmax=1.0
						newmin=0.0
						max=np.max(elbins[1][i])
						min=np.min(elbins[1][i])
						toplot = (newmax-newmin)/(max-min)*(elbins[1][i]-max)+newmax
						toplot = toplot - skydip_fitted
						if i in skydip_mask:
							plt.plot(elbins[0][i],toplot,'b',alpha=0.5)
						else:
							plt.plot(elbins[0][i],toplot,'r',alpha=0.5)

					plt.xlabel('Elevation')
					plt.ylabel('Power (fit subtracted)')
					plt.savefig(self.outdir+'/'+prefix+'/average_combine_sub_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.pdf')
					plt.close()
					plt.clf()


		return

	# Originally by Roger, 'read_sci_FTGI_multi2018.py'
	def get_sci(self, filenames,quiet=False):
		numfiles = len(filenames)
		# data=np.empty(shape=(124,60*4000*numfiles),dtype=float)
		dat=np.empty(4000,dtype='i8')
		data = []
		for i in range(0,124):
			data.append([])
		for j in range(0,numfiles):
			print(filenames[j])
			with open(filenames[j], "rb") as f:
				for k in range(60):
					for i in range(124):
						campo1= np.fromfile(f,count=1, dtype='>a2')
						campo2= np.fromfile(f,count=1, dtype='>u2')            
						campo3= np.fromfile(f,count=1, dtype='>u2')
						campo4= np.fromfile(f,count=1, dtype='>u1')
						campo5= np.fromfile(f,count=1, dtype='>u1')
						campo6= np.fromfile(f,count=1, dtype='>u1')
						campo7= np.fromfile(f,count=1, dtype='>u1')
						Tstamp= np.fromfile(f,count=1, dtype='>u4')
						Time=   np.fromfile(f,count=1, dtype='>u4')
						SID=    np.fromfile(f,count=1, dtype='>u2')
						AqErr=  np.fromfile(f,count=1, dtype='>u1')
						channel_ID=  np.fromfile(f,count=1, dtype='>u1')
						cal_flag=  np.fromfile(f,count=1, dtype='>u1')
						cal_sw=  np.fromfile(f,count=1, dtype='>u1')
						ph_sw=  np.fromfile(f,count=1, dtype='>u1') 
						repuesto1=  np.fromfile(f,count=1, dtype='>u1') 
						Nphstates=  np.fromfile(f,count=1, dtype='>u4')
						PHseq= np.fromfile(f,count=16, dtype='>u1')
						Processtype= np.fromfile(f,count=1, dtype='>u1')
						repuesto2= np.fromfile(f,count=1, dtype='>u1')
						Pack_count=  np.fromfile(f,count=1, dtype='>u2')
						Sc_factor=np.fromfile(f,count=1, dtype='>f8')
						Samprate=  np.fromfile(f,count=1, dtype='>u4')
						NSample=  np.fromfile(f,count=1, dtype='u4')
						dat= np.fromfile(f,count=4000, dtype='>i4')
						dat=dat*Sc_factor
						if len(dat) != 0:
							# data[i,(j*240000)+(k*4000):(j*240000)+(k*4000)+4000]=dat
							data[i].append(dat)
						if not quiet:
							print(j,k,i)
							print('campo1 = ' + str(campo1))
							print('campo2 = ' + str(campo2))
							print('campo3 = ' + str(campo3))
							print('campo4 = ' + str(campo4))
							print('campo5 = ' + str(campo5))
							print('campo6 = ' + str(campo6))
							print('campo7 = ' + str(campo7))
							print('Tstamp = ' + str(Tstamp))
							print('Time = ' + str(Time))
							print('SID = ' + str(SID))
							print('AqErr = ' + str(AqErr)) 
							print('channel_ID = ' + str(channel_ID))
							print('cal_flag = ' + str(cal_flag))
							print('cal__sw = ' + str(cal_sw))
							print('ph__sw = ' + str(ph_sw) + ' (1- 16KHz, 2- 8KHz)')
							print('repuesto1 = ' + str(repuesto1))
							print('Nphstates = ' + str(Nphstates))
							print('Phase sequence = ' + str(PHseq))
							print('Process type = ' + str(Processtype))
							print('repuesto2 = ' + str(repuesto2))
							print('packet counter = ' + str(Pack_count))
							print('Scale factor = ' + str(Sc_factor))
							print('Sampling rate = ' + str(Samprate))
							print('No of samples = ' + str(NSample))
							print(dat[0:10])
		return np.array(data)

	# Originally by Roger, 'plot_cal_data_sci_FGI_multi2018.py'
	def plot_sci(self, data,channels=range(0,30),pixels=range(0,30),fix_neg=False,offset=True):
		dt2=data
		ds=int(data.size/(8*124))
		# fls=int(ds/30000)
		# DAS=int(input("what DAS channel would you like to calibrate?"))
		for channel in range(0,len(channels)):
			DAS = channels[channel]
			pixel = pixels[channel]
			chan=DAS*4-4
			# reshape array into phase sum cycles (8 4ON/4OFF) and then take Cal OFF from CAL ON at 4KHz/125Hz
			# There is an ambiguity in the order of the CALON/CALOFF. It changes sign after 4 s
			dt3=np.reshape(dt2,(124,ds,8))
			if pixel >= 40:
				# FGI
				if offset:
					adj=1/dt3[chan:chan+4,0,0:4]
					dt4=dt3[:,:,4:]/dt3[:,:,0:4]
					dt3null=dt3[:,:,0:4]/dt3[:,:,0:4]
					dt3[:,:,0:4]=dt3null[:,:,:]
					dt3[:,:,4:]=dt4[:,:,:]
					# dt4=dt4[:,:,:]-dt3null[:,:,:]
					dt4=dt3[:,:,0:4]-dt3[:,:,4:]
				else:
					dt4=dt3[:,:,0:4]-dt3[:,:,4:]
			else:
				# TGI
				if offset:
					adj=1/dt3[chan:chan+4,0,4:]
					dt4=dt3[:,:,0:4]/dt3[:,:,4:]
					dt3null=dt3[:,:,4:]/dt3[:,:,4:]
					dt3[:,:,4:]=dt3null[:,:,:]
					dt3[:,:,0:4]=dt4[:,:,:]
					# dt4=dt4[:,:,:]-dt3null[:,:,:]
					dt4=dt3[:,:,4:]-dt3[:,:,0:4]
				else:
					dt4=dt3[:,:,4:]-dt3[:,:,0:4]
			# Calculate the mean values over a given phase state for all phase states in the file
			phsum0=dt4[chan:chan+4,:,0]
			if pixel >= 40:
				# FGI
				phsum90=dt4[chan:chan+4,:,1] 
				phsum180=dt4[chan:chan+4,:,2]  
			else:
				# TGI
				phsum180=dt4[chan:chan+4,:,1] 
				phsum90=dt4[chan:chan+4,:,2]  

			phsum270=dt4[chan:chan+4,:,3]
			# Calculate I, Q and U
			I=(phsum0+phsum180+phsum90+phsum270)/2
			Q=phsum0-phsum180
			U=phsum90-phsum270
			# Calculate the polar magnitude
			L=np.absolute(Q+U*1j)
			# Caculate the polar angle
			# L_ang=np.angle(Q+U*1j,deg=1)
			L_ang=(0.5*np.arctan2(U,Q)* 180 / np.pi)+90.

			# Do some plots of those
			# print(np.shape(I))
			# print(len(I))
			# for j in range(0,4):
			# 	plt.plot(I[j,0:1000])
			# 	plt.plot(Q[j,0:1000])
			# 	plt.plot(U[j,0:1000])
			# 	plt.savefig('plot_IQU_'+str(DAS)+'_'+str(j)+'.png')
			# 	plt.clf()

			# fft_w = scipy.fftpack.rfft(I[0,:])
			# fft_f = scipy.fftpack.rfftfreq(len(I[0,:]), 1)
			# fft_spectrum = fft_w**2
			# plt.plot(fft_f, fft_w)
			# plt.xscale("log")
			# plt.yscale("log")
			# plt.savefig('plot_I_'+str(DAS)+'_fft.png')
			# plt.clf()

			# Auto-create a mask for the different states
			mask = np.ones(len(L_ang[0,:]))
			for i in range(1,len(L_ang[0,:])):
				if np.abs(L_ang[0,i] - L_ang[0,i-100]) > 1.0 or np.abs(L_ang[1,i] - L_ang[1,i-100]) > 1.0:
					# We have a step. Mark it in the mask.
					mask[i] = 0
			# Ignore any jumps in the first 1,000 samples
			mask[0:1000] = 0
			# ... and the last 1,000 samples
			mask[-1000:-1] = 0
			# Broaden the borders
			testlen = 1500
			for i in range(0,len((mask/testlen)-testlen)):
				if np.sum(np.abs(mask[i*testlen:i*testlen+testlen])) != testlen:
					mask[i*testlen:i*testlen+testlen] = 0

			# Count how many different sections we have, and update the mask so we can extract them.
			num_sections = 1
			notzero = 0
			borders = []
			for i in range(0,len(mask)):
				if mask[i] != 0:
					mask[i] = mask[i] * num_sections
					if notzero == 0:
						borders.append(i)
					notzero = 1
				else:
					if notzero == 1:
						num_sections = num_sections+1
						notzero = 0
						borders.append(i)
			# print(num_sections)

			# print(len(mask))
			# print(np.sum(mask))

			# Calculate file statistics in terms of polar magnitude and angle
			# I_mean=np.empty(shape=(4,fls),dtype=float)
			# L_mean=np.empty(shape=(4,fls),dtype=float)
			# L_mean_rms=np.empty(shape=(4,fls),dtype=float)
			# L_ang_mean=np.empty(shape=(4,fls),dtype=float)
			# L_ang_rms=np.empty(shape=(4,fls),dtype=float)
			I2_mean=np.empty(shape=(4,num_sections-1),dtype=float)
			L2_mean=np.empty(shape=(4,num_sections-1),dtype=float)
			L2_mean_rms=np.empty(shape=(4,num_sections-1),dtype=float)
			L2_ang_mean=np.empty(shape=(4,num_sections-1),dtype=float)
			L2_ang_rms=np.empty(shape=(4,num_sections-1),dtype=float)
			# for k in range(fls):
			# 	I_mean[:,k]=np.mean(I[:,(30000*k+10000):(30000*k+20000)],axis=1)
			# 	L_mean[:,k]=np.mean(L[:,(30000*k+10000):(30000*k+20000)],axis=1)
			# 	L_mean_rms[:,k]=np.std(L[:,(30000*k+10000):(30000*k+20000)],axis=1)
			# 	L_ang_mean[:,k]=np.mean(L_ang[:,(30000*k+10000):(30000*k+20000)],axis=1)
			# 	L_ang_rms[:,k]=np.std(L_ang[:,(30000*k+10000):(30000*k+20000)],axis=1)
			# Print out in table format using pandas
			for j in range(0,num_sections-1):
				# stepmask = np.zeros(len(mask))
				# stepmask[mask == j+1] = 1
				# print(np.sum(stepmask))
				# print(np.mean(I[:,mask==j+1],axis=1))
				I2_mean[:,j] = np.mean(I[:,mask==j+1],axis=1)
				L2_mean[:,j] = np.mean(L[:,mask==j+1],axis=1)
				L2_mean_rms[:,j] = np.std(L[:,mask==j+1],axis=1)
				if fix_neg:
					L2_ang_mean[:,j] = np.mean(-L_ang[:,mask==j+1],axis=1)
				else:
					L2_ang_mean[:,j] = np.mean(L_ang[:,mask==j+1],axis=1)
				L2_ang_rms[:,j] = np.std(L_ang[:,mask==j+1],axis=1)

			# L2_ang_mean[:,:] = L2_ang_mean[:,:] + 90.0
			L2_ang_mean[L2_ang_mean > 90.0] = L2_ang_mean[L2_ang_mean > 90.0] - 180.0
			L2_ang_mean[L2_ang_mean < -90.0] = L2_ang_mean[L2_ang_mean < -90.0] + 180.0

			for k in range(4):
				# print(num_sections)
				if num_sections == 6:
					p_ang=["-45º","-22.5º","0º","22.5º","45º"]
					p_ang_num = [-45,-22.5,0,22.5,45]
					v_out=[k+1,k+1,k+1,k+1,k+1]
					maxnum = 5
				else:
					p_ang=["-45º","-22.5º","0º","22.5º","45º","67.5º"]
					p_ang_num = [-45,-22.5,0,22.5,45,67.5]
					v_out=[k+1,k+1,k+1,k+1,k+1,k+1]
					maxnum = 6
				# P_stats=list(zip(v_out,p_ang,L_mean[k,:],L_mean_rms[k,:],L_ang_mean[k,:],L_ang_rms[k,:]))
				# P_stats_table=pd.DataFrame(data=P_stats,columns=["Vout","Output","P mag (volts)","P mag rms","P angle (deg)","P angle rms"])
				# print(P_stats_table)  
				p_ang_corr = L2_ang_mean[k,0:maxnum]+p_ang_num
				p_ang_corr[p_ang_corr > 90.0] = p_ang_corr[p_ang_corr > 90.0] - 180.0
				p_ang_corr[p_ang_corr < -90.0] = p_ang_corr[p_ang_corr < -90.0] + 180.0

				P2_stats=list(zip(v_out,p_ang,L2_mean[k,:],L2_mean_rms[k,:],L2_ang_mean[k,:],L2_ang_rms[k,:],p_ang_corr,np.abs(L2_mean[k,:]*100/I2_mean[k,:])))
				P2_stats_table=pd.DataFrame(data=P2_stats,columns=["Vout","Output","P mag (volts)","P mag rms","P angle (deg)","P angle rms","P angle (deg), corrected","%pol"])
				print(P2_stats_table)
				if offset:
				    print(adj[k,:])
				print("Average P ang (corrected): " + '{:3.1f}'.format(np.mean(p_ang_corr)))
				print("Average %pol: " + '{:3.1f}'.format(np.mean(np.abs(L2_mean[k,:]*100/I2_mean[k,:]))))
			# Plot data on screen and into a pdf file
			fig=plt.figure()
			fig.suptitle('Scientific data for Pixel '+str(pixel))
			plt.subplot(2,2,1)
			tim=(np.arange(4)+1)*2
			tim1=np.arange(8)+1
			tim2=np.arange(ds)/125              
			plt.grid(True)
			if pixel >= 40:
				plt.xticks([0,1,2,3,4,5,6,7],('0','90','180','270','0','90','180','270'))
			else:
				plt.xticks([0,1,2,3,4,5,6,7],('0','180','90','270','0','180','90','270'))
			plt.plot(dt3[chan,0,:],label='V1')
			plt.plot(dt3[chan+1,0,:],label='V2')
			plt.plot(dt3[chan+2,0,:],label='V3')
			plt.plot(dt3[chan+3,0,:],label='V4')
			plt.legend()
			plt.ylabel('voltage output, V')
			plt.xlabel('Phase Switch angle, deg')
			plt.title('CalOn/CalOff Cycle')
			plt.subplot(2,2,2)
			plt.grid(True)
			if pixel > 40:
				plt.xticks([0,1,2,3],('0','90','180','270'))
			else:
				plt.xticks([0,1,2,3],('0','180','90','270'))
			plt.plot(dt4[chan,0,:])
			plt.plot(dt4[chan+1,0,:])
			plt.plot(dt4[chan+2,0,:])
			plt.plot(dt4[chan+3,0,:])
			plt.ylabel('voltage output, V')
			plt.xlabel('Phase Switch Angle, deg')
			plt.title('Cal polar signal')
			plt.subplot(2,2,3)
			plt.grid(True)
			plt.plot(tim2,L[0,:])
			plt.plot(tim2,L[1,:])
			plt.plot(tim2,L[2,:])
			plt.plot(tim2,L[3,:])
			for i in range(0,len(borders)):
				plt.axvline(x=float(borders[i])*8.0/1000.0)
			plt.title('Detector polar magnitude')
			plt.xlabel('time,s')
			plt.ylabel('voltage, V')
			plt.subplot(2,2,4)
			plt.grid(True)
			plt.yticks(np.arange(0,180,30))
			plt.plot(tim2,L_ang[0,:])
			plt.plot(tim2,L_ang[1,:])
			plt.plot(tim2,L_ang[2,:])
			plt.plot(tim2,L_ang[3,:])
			for i in range(0,len(borders)):
				plt.axvline(x=float(borders[i])*8.0/1000.0)
			# plt.plot(tim2,mask*50)
			plt.title('Detector polar Angle')
			plt.xlabel('time,s')
			plt.ylabel('angle, deg')
			plt.subplots_adjust(bottom=0.1,left=0.1, right=0.9, top=0.8,wspace=0.4,hspace=0.4)
			fig.savefig('plot_sci_'+str(pixel)+'.pdf')
			fig.clf()
		return


	# Originally by Roger, 'read_eng_FTGI_2018.py'
	def get_eng(self, filenames,quiet=False):
		data=np.empty(shape=(8,60*40*4000),dtype=float)
		dat=np.empty(4000,dtype='i8')
		with open(filenames, "rb") as f:
			for k in range(60*40):
				for ch in range(8):
					campo1= np.fromfile(f,count=1, dtype='>a2')
					campo2= np.fromfile(f,count=1, dtype='>u2')            
					campo3= np.fromfile(f,count=1, dtype='>u2')
					campo4= np.fromfile(f,count=1, dtype='>u1')
					campo5= np.fromfile(f,count=1, dtype='>u1')
					campo6= np.fromfile(f,count=1, dtype='>u1')
					campo7= np.fromfile(f,count=1, dtype='>u1')
					Tstamp= np.fromfile(f,count=1, dtype='>u4')
					Time=   np.fromfile(f,count=1, dtype='>u4')
					SID=    np.fromfile(f,count=1, dtype='>u2')
					AqErr=  np.fromfile(f,count=1, dtype='>u1')
					channel_ID=  np.fromfile(f,count=1, dtype='>u1')
					cal_flag=  np.fromfile(f,count=1, dtype='>u1')
					cal_sw=  np.fromfile(f,count=1, dtype='>u1')
					ph_sw=  np.fromfile(f,count=1, dtype='>u1') 
					repuesto1=  np.fromfile(f,count=1, dtype='>u1') 
					Nphstates=  np.fromfile(f,count=1, dtype='>u4')
					PHseq= np.fromfile(f,count=16, dtype='>u1')
					Processtype= np.fromfile(f,count=1, dtype='>u1')
					repuesto2= np.fromfile(f,count=1, dtype='>u1')
					Pack_count=  np.fromfile(f,count=1, dtype='>u2')
					Sc_factor=np.fromfile(f,count=1, dtype='>f8')
					Samprate=  np.fromfile(f,count=1, dtype='>u4')
					NSample=  np.fromfile(f,count=1, dtype='u4')
					dat= np.fromfile(f,count=4000, dtype='>i4')
					dat=dat*Sc_factor
					data[ch,(k*4000):(k*4000)+4000]=dat
					"""
					Correct for calibration offset by taking off the first 80 phase states

					for n in range(8):
					   data[n-1,:]=np.roll(data[n-1,:],-80)          
					   
					"""
					if not quiet:
						print(k)      
						print('campo1 = ' + str(campo1))
						print('campo2 = ' + str(campo2))
						print('campo3 = ' + str(campo3))
						print('campo4 = ' + str(campo4))
						print('campo5 = ' + str(campo5))
						print('campo6 = ' + str(campo6))
						print('campo7 = ' + str(campo7))
						print('Tstamp = ' + str(Tstamp))
						print('Time = ' + str(Time))
						print('SID = ' + str(SID))
						print('AqErr = ' + str(AqErr)) 
						print('channel_ID = ' + str(channel_ID))
						print('cal_flag = ' + str(cal_flag))
						print('cal__sw = ' + str(cal_sw))
						print('ph__sw = ' + str(ph_sw) + ' (1- 16KHz, 2- 8KHz)')
						print('repuesto1 = ' + str(repuesto1))
						print('Nphstates = ' + str(Nphstates))
						print('Phase sequence = ' + str(PHseq))
						print('Process type = ' + str(Processtype))
						print('repuesto2 = ' + str(repuesto2))
						print('packet counter = ' + str(Pack_count))
						print('Scale factor = ' + str(Sc_factor))
						print('Sampling rate = ' + str(Samprate))
						print('No of samples = ' + str(NSample))
						print(dat[0:10])
		return data

	# Originally by Roger, 'plot_cal_data_sci_FGI_multi2018.py'
	def plot_eng(self, data,pixel=0):
		dt2=data
		# ns is the positive transition samples and ne the negative slope transition samples
		# that have been removed from the calculation
		ns=1
		ne=1
		# reshape array into phase cycles (640) and then take Cal OFF from CAL ON at 4KHz/125Hz
		dt3=np.reshape(dt2,(8,-1,1280))
		# (DT4 calculation moved below)
		# Reduce the phase states from 16 to 4 by summing equivalent states taking the positive
		# and negative transition samples out of the sumation
		if pixel >= 40:
			# FGI
			dt4=dt3[:,:,0:640]-dt3[:,:,640:]
			ph0=(dt4[:,:,0+ns:40-ne]+dt4[:,:,240+ns:280-ne]+dt4[:,:,400+ns:440-ne]+dt4[:,:,480+ns:520-ne])/4
			ph90=(dt4[:,:,40+ns:80-ne]+dt4[:,:,200+ns:240-ne]+dt4[:,:,440+ns:480-ne]+dt4[:,:,600+ns:640-ne])/4
			ph180=(dt4[:,:,120+ns:160-ne]+dt4[:,:,160+ns:200-ne]+dt4[:,:,320+ns:360-ne]+dt4[:,:,560+ns:600-ne])/4 
			ph270=(dt4[:,:,80+ns:120-ne]+dt4[:,:,280+ns:320-ne]+dt4[:,:,360+ns:400-ne]+dt4[:,:,520+ns:560-ne])/4
		else:
			# TGI
			dt4=dt3[:,:,640:]-dt3[:,:,0:640]
			ph0=(dt4[:,:,0+ns:40-ne]+dt4[:,:,240+ns:280-ne]+dt4[:,:,400+ns:440-ne]+dt4[:,:,480+ns:520-ne])/4
			ph180=(dt4[:,:,40+ns:80-ne]+dt4[:,:,280+ns:320-ne]+dt4[:,:,440+ns:480-ne]+dt4[:,:,520+ns:560-ne])/4
			ph90=(dt4[:,:,120+ns:160-ne]+dt4[:,:,200+ns:240-ne]+dt4[:,:,320+ns:360-ne]+dt4[:,:,560+ns:600-ne])/4 
			ph270=(dt4[:,:,80+ns:120-ne]+dt4[:,:,160+ns:200-ne]+dt4[:,:,360+ns:400-ne]+dt4[:,:,600+ns:640-ne])/4

		# Calculate the mean values over a given phase state for all phase states in the file (equal to scientific data)
		phsum0=np.mean(ph0,axis=2)
		phsum90=np.mean(ph90,axis=2) 
		phsum180=np.mean(ph180,axis=2)  
		phsum270=np.mean(ph270,axis=2)
		# Calculate Q and U
		I=(phsum0+phsum180+phsum90+phsum270)/2
		Q=phsum0-phsum180
		U=phsum90-phsum270
		# Calculate the polar magnitude
		L=np.absolute(Q+U*1j)
		# Calculate the polar angle
		# L_ang=np.angle(Q+U*1j,deg=1)
		L_ang=(0.5*np.arctan2(U,Q)* 180 / np.pi)+90.
		# Calculate file statistics in terms of polar magnitude and angle
		I_mean=np.mean(I[0:4,:],axis=1)
		L_mean=np.mean(L[0:4,:],axis=1)
		L_mean_rms=np.std(L[0:4,:],axis=1)
		L_ang_mean=np.mean(L_ang[0:4,500:6500],axis=1)
		# If over 90 degrees, subtract 180 from it.
		L_ang_mean[L_ang_mean > 90.0] = L_ang_mean[L_ang_mean > 90.0] - 180.0

		L_ang_rms=np.std(L_ang[0:4,500:6500],axis=1)
		# Print out in table format using pandas
		V_out=["Vout1","Vout2","Vout3","Vout4"]
		P_stats=list(zip(V_out,L_mean,L_mean_rms,L_ang_mean,L_ang_rms,np.abs(L_mean*100/I_mean)))
		P_stats_table=pd.DataFrame(data=P_stats,columns=["Output","P mag (volts)","P mag rms","P angle (deg)","P angle rms","%pol"])
		print(P_stats_table)
		# Plot data on screen and into a pdf file
		fig=plt.figure()
		fig.suptitle('Engineering data for Pixel ' + str(pixel))
		plt.subplot(2,2,1)
		tim1=np.arange(1280)/160
		tim2=np.arange(7500)/125              
		plt.grid(True)
		plt.plot(tim1,dt3[0,0,:],label='V1')
		plt.plot(tim1,dt3[1,0,:],label='V2')
		plt.plot(tim1,dt3[2,0,:],label='V3')
		plt.plot(tim1,dt3[3,0,:],label='V4')
		plt.legend()
		plt.ylabel('voltage output, V')
		plt.xlabel('time,msecs')
		plt.title('polar signal')
		plt.subplot(2,2,2)
		plt.grid(True)
		if pixel >= 40:
			# FGI
			plt.xticks([0,1,2,3],('0','90','180','270'))
			plt.plot([phsum0[0,0],phsum90[0,0],phsum180[0,0],phsum270[0,0]])
			plt.plot([phsum0[1,0],phsum90[1,0],phsum180[1,0],phsum270[1,0]])
			plt.plot([phsum0[2,0],phsum90[2,0],phsum180[2,0],phsum270[2,0]])
			plt.plot([phsum0[3,0],phsum90[3,0],phsum180[3,0],phsum270[3,0]])
		else:
			# TGI
			plt.xticks([0,1,2,3],('0','180','90','270'))
			plt.plot([phsum0[0,0],phsum180[0,0],phsum90[0,0],phsum270[0,0]])
			plt.plot([phsum0[1,0],phsum180[1,0],phsum90[1,0],phsum270[1,0]])
			plt.plot([phsum0[2,0],phsum180[2,0],phsum90[2,0],phsum270[2,0]])
			plt.plot([phsum0[3,0],phsum180[3,0],phsum90[3,0],phsum270[3,0]])
		plt.ylabel('voltage output, V')
		plt.xlabel('phase state, deg.')
		plt.title('scientific data polar signal')
		plt.subplot(2,2,3)
		plt.grid(True)
		plt.plot(tim2,L[0,:])
		plt.plot(tim2,L[1,:])
		plt.plot(tim2,L[2,:])
		plt.plot(tim2,L[3,:])
		plt.title('Detector polar magnitude')
		plt.xlabel('time,s')
		plt.ylabel('voltage, volts')
		plt.subplot(2,2,4)
		plt.grid(True)
		plt.title('Detector polar Angle')
		plt.plot(tim2,L_ang[0,:])
		plt.plot(tim2,L_ang[1,:])
		plt.plot(tim2,L_ang[2,:])
		plt.plot(tim2,L_ang[3,:])
		plt.xlabel('time,s')
		plt.ylabel('angle, deg.')
		plt.subplots_adjust(bottom=0.1,left=0.1, right=0.9, top=0.8,wspace=0.4,hspace=0.4)
		fig.savefig('plot_eng_'+str(pixel)+'.pdf')
		return
