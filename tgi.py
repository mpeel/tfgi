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

def linfit(x, A, B):
	return A*x+B

def fit_kneefreq(freq, param):
	sigma, fknee, alpha = param
	return sigma**2 * (1 + (fknee / freq)**alpha)

def compute_residuals(param, data, freq):
	model = fit_kneefreq(freq, param)
	residual = np.log(data / model)
	return residual


class tgi:
	def __init__(self,indir="/net/nas/proyectos/quijote2",outdir='',pixelfileloc="/net/nas/proyectos/quijote2/etc/qt2_pixel_masterfile.",pixelposfileloc="/net/nas/proyectos/quijote2/etc/tgi_fgi_horn_positions_table.txt"):
		self.numpixels = 31
		self.indir = indir
		self.outdir = outdir
		self.telescope = EarthLocation(lat=28.300224*u.deg, lon=-16.510113*u.deg, height=2390*u.m)
		self.jd_ref = 2456244.5 # Taken as the first day of the Commissioning (13/Nov/2012, 0.00h)
		self.nside = 512
		self.npix = hp.nside2npix(self.nside)
		self.pixeljds, self.pixelarrays = self.read_pixel_masterfiles(pixelfileloc)
		self.pixelpositions = self.read_pixel_positions(pixelposfileloc)
		self.apobserver = astroplan.Observer(latitude=28.300224*u.deg, longitude=-16.510113*u.deg, elevation=2390*u.m)
		# NB: Galactic doesn't work with hour angles
		self.coordsys = 1 # 0 = Galactic, 1 = ICRS.

	def ensure_dir(self,f):
		os.makedirs(f, exist_ok=True)

	# Read in a single pixel masterfile, and return the dictionary of the results
	def read_pixel_masterfile(self, filename):
		jd = 0
		array = []
		with open(filename) as f:
			for line in f:
				if jd == 0:
					jd = float(line.strip())
				else:
					val = line.strip().split()
					# Convert to integers
					val = [int(x) for x in val]
					array.append(val)
		return jd, array
		
	# Read in a single pixel masterfile, and return the dictionary of the results
	def read_pixel_positions(self, filename):
		array = []
		with open(filename) as f:
			for line in f:
				if '*' not in line:
					val = line.strip().split()
					array.append([int(val[0]), float(val[1]), float(val[2])])
		return array
		
	def read_pixel_masterfiles(self, prefix):
		jds = []
		arrays = []
		for i in range(1,500):
			print(prefix+format(i, '03d')+'.txt')
			try:
				jd, array = self.read_pixel_masterfile(prefix+format(i, '03d')+'.txt')
			except:
				break
			jds.append(jd)
			arrays.append(array)
		# print(jds)
		# exit()
		return jds, arrays

	def get_pixel_info(self, jd, das):
		print(jd)
		print(das)
		print(self.pixeljds)
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
				print(pixelposition)
				return {'fp':line[0],'pixel':line[1],'fem':line[2],'bem':line[3],'das':line[5],'tgi':tgi,'fgi':fgi,'x_pos':pixelposition[1],'y_pos':pixelposition[2]}
				# return line
		return []

	def read_tod(self, prefix, numfiles=50,quiet=True):
		# Read in the data
		jd = np.empty(0)
		az = np.empty(0)
		el = np.empty(0)
		data = np.empty((0,0,0))
		for i in range(0,numfiles):
			print(self.indir+'/'+prefix+'-'+format(i, '04d')+'.tod2')
			try:
				inputfits = fits.open(self.indir+'/'+prefix+'-'+format(i, '04d')+'.tod2')
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
		print(' New data length: ' + str(len(data)/(self.numpixels*4*4)))
		data = data.reshape(4, self.numpixels, 4, ndata, order='F')
		print(' New data shape: ' + str(np.shape(data)))
		return az, el, jd, data

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
			healpix_pixel = hp.ang2pix(self.nside, (3.1415/2)-Angle(skypos.b).radian, Angle(skypos.l).radian)
			pos = (np.median(Angle(skypos.l).degree),np.median((Angle(skypos.b).degree)))
		else:
			print(skypos)
			healpix_pixel = hp.ang2pix(self.nside, (3.1415/2)-Angle(skypos.dec).radian, Angle(skypos.ra).radian)
			pos = (np.median(Angle(skypos.ra).degree),np.median((Angle(skypos.dec).degree)))
		return healpix_pixel, pos

	def plot_tod(self, data, outputname):
		plt.plot(data,'b.')
		plt.xlabel('Samples')
		plt.ylabel('Power')
		plt.savefig(outputname)
		plt.close()
		plt.clf()
		return

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

	def converttopol(self, data, ordering=[0,1,2,3],pa=[]):
		newdata = np.zeros(np.shape(data))
		newdata[0] = (data[ordering[0],:] + data[ordering[1],:])# / 2.0
		newdata[1] = (data[ordering[0],:] - data[ordering[1],:])# / 2.0
		newdata[2] = (data[ordering[2],:] - data[ordering[3],:])# / 2.0
		newdata[3] = (data[ordering[2],:] + data[ordering[3],:])# / 2.0

		if len(pa) > 0:
			# ang = pa
			ang = 2.0*pa #* 3.1415/180.0
			# print(ang)
			Q = newdata[1,:]*np.cos(ang) + newdata[2,:]*np.sin(ang)
			U = -newdata[1,:]*np.sin(ang) + newdata[2,:]*np.cos(ang) 
			newdata[1] = Q
			newdata[2] = U

		return newdata


	def analyse_tod(self, prefix,pixelrange=range(0,31),detrange=range(0,4),phaserange=range(0,4),plotlimit=0.0,quiet=False,dofft=False,plottods=True,dopol=False,numfiles=50):
		print(self.outdir+'/'+prefix)
		self.ensure_dir(self.outdir+'/'+prefix)

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

					# if pix == 3 or pix == 4 or pix == 5 or pix == 21 or pix == 22 or pix == 23 or pix == 24 or pix == 25:
					# print(pix)
					# Get the maximum value in the map
					maxval = np.max(skymap)
					# Get a sample of data from the start of the measurement
					std = np.std(data[det][pix][j][0:4000])
					print(maxval)
					print(std)
					print(std*(400/maxval)/np.sqrt(1000.0))

					hp.write_map(self.outdir+'/'+prefix+'/skymap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.fits',skymap,overwrite=True)
					hp.mollview(skymap)
					plt.savefig(self.outdir+'/'+prefix+'/skymap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
					plt.close()
					plt.clf()
					hp.write_map(self.outdir+'/'+prefix+'/hitmap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.fits',hitmap,overwrite=True)
					hp.mollview(hitmap)
					plt.savefig(self.outdir+'/'+prefix+'/hitmap_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
					hp.gnomview(skymap,rot=centralpos,reso=5)
					plt.savefig(self.outdir+'/'+prefix+'/skymap_zoom_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
					plt.close()
					plt.clf()
					if plotlimit != 0.0:
						hp.mollview(skymap,min=-plotlimit,max=plotlimit)
						plt.savefig(self.outdir+'/'+prefix+'/skymap_cut_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						plt.close()
						plt.clf()
						hp.gnomview(skymap,rot=centralpos,reso=5,min=-plotlimit,max=plotlimit)
						plt.savefig(self.outdir+'/'+prefix+'/skymap_zoom_cut_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
						plt.close()
						plt.clf()


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
		self.ensure_dir(self.outdir+'/'+prefix)

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
			inputmap = hp.read_map(skymaps[i])
			inputhitmap = hp.read_map(hitmaps[i])
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

	# Originally by Roger, 'read_sci_FTGI_2018.py'
	def get_sci(self, filename,quiet=False):
		data=np.empty(shape=(124,60*4000),dtype=float)
		dat=np.empty(4000,dtype='i8')
		with open(filename, "rb") as f:
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
					data[i,(k*4000):(k*4000)+4000]=dat
					
					if not quiet:
						print(k,i)      
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
	def plot_sci(self, data):
		dt2=data
		ds=int(data.size/(8*124))
		fls=int(ds/30000)
		# DAS=int(input("what DAS channel would you like to calibrate?"))
		for DAS in range(1,30):
			chan=DAS*4-4
			"reshape array into phase sum cycles (8 4ON/4OFF) and then take Cal OFF from CAL ON at 4KHz/125Hz"
			"There is an ambiguity in the order of the CALON/CALOFF. It changes sign after 4 s"
			dt3=np.reshape(dt2,(124,ds,8))
			dt4=dt3[:,:,4:]-dt3[:,:,0:4]
			""
			"Calculate the mean values over a given phase state for all phase states in the file"      
			phsum0=dt4[chan:chan+4,:,0]
			phsum180=dt4[chan:chan+4,:,1] 
			phsum90=dt4[chan:chan+4,:,2]  
			phsum270=dt4[chan:chan+4,:,3]
			"Calculate I, Q and U "
			I=(phsum0+phsum180+phsum90+phsum270)/2
			Q=phsum0-phsum180
			U=phsum90-phsum270
			"Calculate the polar magnitude"
			L=np.absolute(Q+U*1j)
			"Caculate the polar angle"
			"L_ang=np.angle(Q+U*1j,deg=1)"
			L_ang=(0.5*np.arctan2(U,Q)* 180 / np.pi)+90.

			# Do some plots of those
			# print(np.shape(I))
			# print(len(I))
			for j in range(0,4):
				plt.plot(I[j,0:1000])
				plt.plot(Q[j,0:1000])
				plt.plot(U[j,0:1000])
				plt.savefig('plot_IQU_'+str(DAS)+'_'+str(j)+'.png')
				plt.clf()

			fft_w = scipy.fftpack.rfft(I[0,:])
			fft_f = scipy.fftpack.rfftfreq(len(I[0,:]), 1)
			fft_spectrum = fft_w**2
			plt.plot(fft_f, fft_w)
			plt.xscale("log")
			plt.yscale("log")
			plt.savefig('plot_I_'+str(DAS)+'_fft.png')
			plt.clf()



			# "Calculate file statistics in terms of polar magnitude and angle"
			# I_mean=np.empty(shape=(4,fls),dtype=float)
			# L_mean=np.empty(shape=(4,fls),dtype=float)
			# L_mean_rms=np.empty(shape=(4,fls),dtype=float)
			# L_ang_mean=np.empty(shape=(4,fls),dtype=float)
			# L_ang_rms=np.empty(shape=(4,fls),dtype=float)
			# for k in range(fls):
			#     I_mean[:,k]=np.mean(I[:,(30000*k+10000):(30000*k+20000)],axis=1)
			#     L_mean[:,k]=np.mean(L[:,(30000*k+10000):(30000*k+20000)],axis=1)
			#     L_mean_rms[:,k]=np.std(L[:,(30000*k+10000):(30000*k+20000)],axis=1)
			#     L_ang_mean[:,k]=np.mean(L_ang[:,(30000*k+10000):(30000*k+20000)],axis=1)
			#     L_ang_rms[:,k]=np.std(L_ang[:,(30000*k+10000):(30000*k+20000)],axis=1)
			# "Print out in table format using pandas"
			# for k in range(4):
			#     p_ang=["-45º","-22.5º","0º","22.5º","45º"]
			#     v_out=[k+1,k+1,k+1,k+1,k+1]
			#     P_stats=list(zip(v_out,p_ang,L_mean[k,:],L_mean_rms[k,:],L_ang_mean[k,:],L_ang_rms[k,:]))
			#     P_stats_table=pd.DataFrame(data=P_stats,columns=["Vout","Output","P mag (volts)","P mag rms","P angle (deg)","P angle rms"])
			#     print(P_stats_table)  
			# "Plot data on screen and into a pdf file"
			# fig=plt.figure()
			# fig.suptitle('Scientific data  for Pixel 24')
			# plt.subplot(2,2,1)
			# tim=(np.arange(4)+1)*2
			# tim1=np.arange(8)+1
			# tim2=np.arange(ds)/125              
			# plt.grid(True)
			# plt.xticks([0,1,2,3,4,5,6,7],('0','180','90','270','0','180','90','270'))
			# plt.plot(dt3[chan,0,:],label='V1')
			# plt.plot(dt3[chan+1,0,:],label='V2')
			# plt.plot(dt3[chan+2,0,:],label='V3')
			# plt.plot(dt3[chan+3,0,:],label='V4')
			# plt.legend()
			# plt.ylabel('voltage output, V')
			# plt.xlabel('Phase Switch angle, deg')
			# plt.title('CalOn/CalOff Cycle')
			# plt.subplot(2,2,2)
			# plt.grid(True)
			# plt.xticks([0,1,2,3],('0','180','90','270'))
			# plt.plot(dt4[chan,0,:])
			# plt.plot(dt4[chan+1,0,:])
			# plt.plot(dt4[chan+2,0,:])
			# plt.plot(dt4[chan+3,0,:])
			# plt.ylabel('voltage output, V')
			# plt.xlabel('Phase Switch Angle, deg')
			# plt.title('Cal polar signal')
			# plt.subplot(2,2,3)
			# plt.grid(True)
			# plt.plot(tim2,L[0,:])
			# plt.plot(tim2,L[1,:])
			# plt.plot(tim2,L[2,:])
			# plt.plot(tim2,L[3,:])
			# plt.title('Detector polar magnitude')
			# plt.xlabel('time,s')
			# plt.ylabel('voltage, V')
			# plt.subplot(2,2,4)
			# plt.grid(True)
			# plt.yticks(np.arange(0,180,30))
			# plt.plot(tim2,L_ang[0,:])
			# plt.plot(tim2,L_ang[1,:])
			# plt.plot(tim2,L_ang[2,:])
			# plt.plot(tim2,L_ang[3,:])
			# plt.title('Detector polar Angle')
			# plt.xlabel('time,s')
			# plt.ylabel('angle, deg')
			# plt.subplots_adjust(bottom=0.1,left=0.1, right=0.9, top=0.8,wspace=0.4,hspace=0.4)
			# fig.savefig('plot_'+str(DAS)+'.pdf')
			# fig.clf()
		return
