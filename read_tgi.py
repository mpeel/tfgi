#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Read in raw TGI data
# 
# Version history:
#
# 02-Apr-2019  M. Peel       Started, based on 
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import pandas as pd
import scipy.fftpack
import matplotlib as mpl
from scipy import signal, optimize
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
from astropy.time import Time
from astropy.coordinates import Angle
from scipy.optimize import curve_fit

mpl.rcParams['agg.path.chunksize'] = 10000

def linfit(x, A, B):
	return A*x+B

def fit_kneefreq(freq, param):
	sigma, fknee, alpha = param
	return sigma**2 * (1 + (fknee / freq)**alpha)

def compute_residuals(param, data, freq):
	model = fit_kneefreq(freq, param)
	residual = np.log(data / model)
	return residual

def get_sci(filename):
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

def plot_sci(data):
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

def analyse_tod(prefix, numfiles=1,pixelrange=range(0,31),detrange=range(0,4),phaserange=range(0,4)):
	# Read in the data
	jd = np.empty(0)
	az = np.empty(0)
	el = np.empty(0)
	data = np.empty((0,0,0))
	for i in range(0,numfiles):
		inputfits = fits.open(prefix+'-'+format(i, '04d')+'.tod2')
		#print(inputfits.info())
		#print(len(inputfits))
		cols = inputfits[1].columns
		col_names = cols.names
		#print(col_names)
		ndata = len(inputfits[1].data.field(0)[0][:])
		nsamples = ndata//4
		if nsamples*4 != ndata:
			print('Oddity in ' + prefix+'-'+format(i, '04d')+'.tod2' + ' - ndata = ' + str(ndata) + ' is not dividable by 4')
		#print(ndata)
		if i == 0:
			jd = inputfits[1].data.field(0)[0][:nsamples*4]
			az = inputfits[1].data.field(1)[0][:nsamples*4]
			el = inputfits[1].data.field(2)[0][:nsamples*4]
		else:
			jd = np.append(jd, inputfits[1].data.field(0)[0][:nsamples*4])
			az = np.append(az,inputfits[1].data.field(1)[0][:nsamples*4])
			el = np.append(el,inputfits[1].data.field(2)[0][:nsamples*4])
		#print(jd[0:100])
		print('Start time: ' + str(np.min(inputfits[1].data.field(0)[0][:nsamples*4])))
		print('End time: ' + str(np.max(inputfits[1].data.field(0)[0][:nsamples*4])))
		print('Duration: ' + str((np.max(inputfits[1].data.field(0)[0][:nsamples*4])-np.min(inputfits[1].data.field(0)[0][:nsamples*4]))*24*60*60) + ' seconds')
		print('There are ' + str(ndata) + " datapoints")
		print('There are ' + str(len(jd)) + " datapoints")	
		print('There are ' + str(len(az)) + " datapoints")	
		#print(az[0:100])
		print('Az range: ' + str(np.min(az)) + ' to ' + str(np.max(az)))
		print('El range: ' + str(np.min(el)) + ' to ' + str(np.max(el)))
		rawdata = inputfits[1].data.field(3)[0][:nsamples*4*4*31]
		print('Raw data array is ' + str(np.shape(rawdata)))
		data = np.append(data, rawdata)
		print(np.shape(data))

	# Reshape the data
	ndata = len(jd)//4
	az = az.reshape(4,ndata,order='F')
	el = el.reshape(4,ndata,order='F')
	jd = jd.reshape(4,ndata,order='F')
	print('JD shape: ' + str(np.shape(jd)))
	print('Az shape: ' + str(np.shape(az)))
	print('Data shape: ' + str(np.shape(data)))
	print(len(data)/(31*4*4))
	numpixels = 31
	data = data.reshape(4, numpixels,4,ndata,order='F')
	print('Data shape: ' + str(np.shape(data)))

	# Calculate the healpix pixel numbers
	telescope = EarthLocation(lat=28.300224*u.deg, lon=-16.510113*u.deg, height=2390*u.m)
	jd_ref = 2456244.5 # Taken as the first day of the Commissioning (13/Nov/2012, 0.00h)

	nside = 512
	npix = hp.nside2npix(nside)
	print(npix)

	# Convert the az/el in the data to healpix pixel numbers
	time = Time(jd[0]+jd_ref, format='jd')
	position = AltAz(az=az[0]*u.deg,alt=el[0]*u.deg,location=telescope,obstime=time)
	skypos = position.transform_to(Galactic)
	healpix_pixel = hp.ang2pix(nside, (3.1415/2)-Angle(skypos.b).radian, Angle(skypos.l).radian)

	# Make maps for each pixel
	for pix in pixelrange:
		maxdata = 50000
		for det in detrange:
			print('Pixel ' + str(pix+1) + ', detector ' + str(det+1))
			plt.plot(data[det][pix][0][:],'b.')
			plt.savefig('plots/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_pre.png')
			plt.clf()
			plt.plot(data[det][pix][0][0:5000],'b.')
			plt.savefig('plots/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_pre_zoom.png')
			plt.clf()

			for j in phaserange:

				fft_w = np.fft.rfft(data[det][pix][j][0:maxdata])
				fft_f = np.fft.rfftfreq(len(data[det][pix][j][0:maxdata]), 1/4000)
				fft_spectrum = np.abs(fft_w)
				#fft_spectrum = np.abs(fft_w)**2
				plt.plot(fft_f, fft_spectrum)
				plt.xscale("log")
				plt.yscale("log")
				plt.savefig('plots/plot_fft_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'.png')
				plt.clf()
				plt.plot(fft_f, fft_spectrum)
				plt.yscale("log")
				plt.savefig('plots/plot_fft_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_lin.png')
				plt.clf()

				# Do a smoothed version
				numbins = 50
				bins = np.linspace(np.log(min(fft_f[1:-1])),np.log(max(fft_f[1:-1])),numbins,endpoint=False)
				values = np.zeros(numbins)
				values2 = np.zeros(numbins)
				matches = np.digitize(np.log(fft_f),bins)
				binmask = np.ones(numbins)
				for i in range(0,numbins):
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
				plt.savefig('plots/plot_fft_'+str(pix+1)+'_'+str(det+1)+'_'+str(j+1)+'_fit.png')
				plt.close()
				plt.clf()

			# Crude baseline removal
			navg = 1500
			option = 0
			for i in range(0,npix//navg):
				start = i*navg
				end = ((i+1)*navg)
				# print('start:' + str(start) + ', end: ' + str(end))
				if option == 0:
					med = np.median(data[det][pix][0][start:end])
					data[det][pix][0][start:end] = data[det][pix][0][start:end] - med
				else:
					xline = range(0,len(data[det][pix][0][start:end]))
					if len(xline) != 0:
						A,B=curve_fit(linfit,xline,data[det][pix][0][start:end])[0]
						# print(A,B)
						data[det][pix][0][start:end] = data[det][pix][0][start:end] - (A*xline+B)

			if option == 0:
				data[det][pix][0][(npix//navg)*navg-1:] = data[det][pix][0][(npix//navg)*navg-1:] - np.median(data[det][pix][0][(npix//navg)*navg-1:])
			else:
				xline = range(0,len(data[det][pix][0][(npix//navg)*navg-1:]))
				if len(xline) != 0:
					A,B=curve_fit(linfit,xline,data[det][pix][0][(npix//navg)*navg-1:])[0]
					# print(A,B)
					data[det][pix][0][(npix//navg)*navg-1:] = data[det][pix][0][(npix//navg)*navg-1:] - (A*xline+B)

			plt.plot(data[det][pix][0][:],'b.')
			plt.savefig('plots/plot_tod_'+str(pix+1)+'_'+str(det+1)+'.png')
			plt.clf()
			plt.plot(data[det][pix][0][0:5000],'b.')
			plt.savefig('plots/plot_tod_'+str(pix+1)+'_'+str(det+1)+'_zoom.png')
			plt.clf()

			skymap = np.zeros(npix, dtype=np.float)
			hitmap = np.zeros(npix, dtype=np.float)
			for i in range(0,len(healpix_pixel)):
				skymap[healpix_pixel[i]] = skymap[healpix_pixel[i]] + data[det][pix][0][i]
				hitmap[healpix_pixel[i]] = hitmap[healpix_pixel[i]] + 1
			for i in range(0,len(skymap)):
				if hitmap[i] >= 1:
					skymap[i] = skymap[i]/hitmap[i]
				else:
					skymap[i] = hp.pixelfunc.UNSEEN

			if pix == 3 or pix == 4 or pix == 5 or pix == 21 or pix == 22 or pix == 23 or pix == 24 or pix == 25:
				# print(pix)
				# Get the maximum value in the map
				maxval = np.max(skymap)
				# Get a sample of data from the start of the measurement
				std = np.std(data[det][pix][0][0:4000])
				print(maxval)
				print(std)
				print(std*(400/maxval)/np.sqrt(1000.0))

			hp.write_map('maps/skymap_'+str(pix+1)+'_'+str(det+1)+'.fits',skymap,overwrite=True)
			hp.mollview(skymap)
			plt.savefig('plots/skymap_'+str(pix+1)+'_'+str(det+1)+'.png')
			plt.clf()
			hp.write_map('maps/hitmap_'+str(pix+1)+'_'+str(det+1)+'.fits',hitmap,overwrite=True)
			hp.mollview(hitmap)
			plt.savefig('plots/hitmap_'+str(pix+1)+'_'+str(det+1)+'.png')
			hp.gnomview(skymap,rot=(185,-7),reso=4)
			plt.savefig('plots/skymap_crab_'+str(pix+1)+'_'+str(det+1)+'.png')
			plt.clf()
	exit()

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


# data = get_sci('testdata/CRAB-19-03-01-18-07-38-0000.sci2')
# plot_sci(data)

# print(data)
# exit()
# dt = np.dtype([('time', [('min', int), ('sec', int)]),('temp', float)])

#11 March 2019, 1728

#analyse_tod('testdata/CRAB-190311-1728',numfiles=15)
#analyse_tod('testdata/CRAB-190409-2029',numfiles=14)
# analyse_tod('testdata/CRAB-190410-1529',numfiles=14)
analyse_tod('testdata/CRAB-190410-2025',numfiles=14,pixelrange=[3,4,5,21,22,23,24])

