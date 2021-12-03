def read_tod_files(indir, prefix, numpixels, numfiles=150,quiet=True):
	# We need to do two passes - one to get the data sizes, one to actually read in the data
	tot_files = 0
	tot_ndata = 0
	tot_nsamples = 0
	numfiles = 10
	for i in range(0,numfiles):
		print(indir+'/'+prefix+'-'+format(i, '04d')+'.tod2')
		try:
			inputfits = fits.open(indir+'/'+prefix+'-'+format(i, '04d')+'.tod2')
		except:
			break
		tot_files += 1
		ndata = len(inputfits[1].data.field(0)[0][:])
		tot_ndata += ndata
		tot_nsamples += ndata//4
		inputfits.close()
	print('Total: ' + str(tot_files) + ' files, ' + str(tot_ndata) + ' data, ' + str(tot_nsamples) + ' samples.')

	# Read in the data
	jd = np.empty(tot_ndata)
	az = np.empty(tot_ndata)
	el = np.empty(tot_ndata)
	data = np.empty((4, numpixels, 4, tot_ndata))
	start_index = 0
	for i in range(0,tot_files):
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
			ndata = nsamples*4
		jd[start_index:start_index+ndata] = inputfits[1].data.field(0)[0][:ndata]
		az[start_index:start_index+ndata] = inputfits[1].data.field(1)[0][:ndata]
		el[start_index:start_index+ndata] = inputfits[1].data.field(2)[0][:ndata]
		rawdata = inputfits[1].data.field(3)[0][:ndata*4*31]
		rawdata = rawdata[0:ndata].T
		# print(ndata)
		# print(np.shape(rawdata))
		# print(np.shape(rawdata[0:ndata]))
		# rawdata = rawdata.reshape(ndata*124,order='C')
		# rawdata = rawdata.reshape(4, numpixels, 4, ndata, order='F')
		# print(nsamples*4*4*31)
		# rawdata = rawdata[:nsamples*4*4*31]
		for j in range(0,4):
			for k in range(0,numpixels):
				for l in range(0,4):
					print(j*k+l)
					data[j][k][l][start_index:start_index+ndata] = rawdata[j*k+l][0:ndata]

		if not quiet:
			print(' Start time: ' + str(np.min(inputfits[1].data.field(0)[0][:nsamples*4])))
			print(' End time: ' + str(np.max(inputfits[1].data.field(0)[0][:nsamples*4])))
			print(' Duration: ' + str((np.max(inputfits[1].data.field(0)[0][:nsamples*4])-np.min(inputfits[1].data.field(0)[0][:nsamples*4]))*24*60*60) + ' seconds')
			print(' There are ' + str(ndata) + " datapoints")
			print(' Az range: ' + str(np.min(inputfits[1].data.field(1)[0][:nsamples*4])) + ' to ' + str(np.max(inputfits[1].data.field(1)[0][:nsamples*4])))
			print(' El range: ' + str(np.min(inputfits[1].data.field(2)[0][:nsamples*4])) + ' to ' + str(np.max(inputfits[1].data.field(2)[0][:nsamples*4])))
			print(' Raw data array is ' + str(np.shape(rawdata)))
		start_index += ndata
		inputfits.close()

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
	# print(' New data length: ' + str(len(data)/(numpixels*4*4)))
	# data = data.reshape(4, numpixels, 4, ndata, order='F')
	# print(' New data shape: ' + str(np.shape(data)))
	return az, el, jd, data
