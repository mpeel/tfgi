import numpy as np
import matplotlib.pyplot as plt

files = ['MOON-190411-1525','CRAB-190410-1529','CRAB-190410-2025','MOON-190223-0528','MOON-190224-0528','MOON-190225-0328']
for file in files:
	array = np.loadtxt('output/'+file+'/_gauflux.txt')
	x = range(0,len(array[:,4]))
	plt.plot(x,array[:,4],label=file)

l = plt.legend(prop={'size':6})
l.set_zorder(20)
plt.savefig('gauflux.pdf')

plt.clf()
conv_from = [3,4,5,21,22,23,24,25,26]

for file in files:
	array = np.loadtxt('output/'+file+'/_gauflux.txt')
	calval = 0.0
	print(len(array))
	# todelete = np.where(array[:,0] == 23)
	# array = np.delete(array, todelete, axis=0)
	# todelete = np.where(array[:,0] == 25)
	# array = np.delete(array, todelete, axis=0)
	# todelete = np.where(array[:,0] == 26)
	# array = np.delete(array, todelete, axis=0)
	for i in range(0,len(array[:,4])):
		try:
			# if array[i,0] == 23 or array[i,0] == 26:
			# 	array = np.delete(array, (i), axis=0)
			# else:
			array[i,0] = conv_from.index(array[i,0]-1)*16 + 4*(array[i,1]-1) + (array[i,2]-1)
			if array[i,2] == 1:
				calval = array[i,4]
			array[i,4] = array[i,4] / calval
		except:
			continue
	# print(array[:,4])
	print(len(array))
	# x = range(0,len(array[:,4]))
	plt.plot(array[:,0],array[:,4],label=file)

l = plt.legend(prop={'size':6})
l.set_zorder(20)
plt.savefig('gauflux_r.pdf')
plt.clf()

for file in files:
	array = np.loadtxt('output/'+file+'/_gauflux.txt')
	calval = 0.0
	print(len(array))
	todelete = np.where(array[:,0] == 23)
	array = np.delete(array, todelete, axis=0)
	todelete = np.where(array[:,0] == 25)
	array = np.delete(array, todelete, axis=0)
	todelete = np.where(array[:,0] == 26)
	array = np.delete(array, todelete, axis=0)
	array[:,4] = array[:,4] * array[:,5] * array[:,6]
	for i in range(0,len(array[:,4])):
		try:
			# if array[i,0] == 23 or array[i,0] == 26:
			# 	array = np.delete(array, (i), axis=0)
			# else:
			if array[i,2] == 1:
				calval = array[i,4]*array[i,5]*array[i,6]
			array[i,4] = array[i,4]*array[i,5]*array[i,6] / calval
		except:
			continue
	# print(array[:,4])
	print(len(array))
	x = range(0,len(array[:,4]))
	plt.plot(x,array[:,4],label=file)

plt.ylim((0.98,1.02))
l = plt.legend(prop={'size':6})
l.set_zorder(20)
plt.savefig('gauflux_r_scale.pdf')


files = ['MOON-190411-1525','CRAB-190410-1529','CRAB-190410-2025','MOON-190223-0528','MOON-190224-0528','MOON-190225-0328']
for file in files:
	array = np.loadtxt('output/'+file+'/_gauflux.txt')
	x = range(0,len(array[:,4]))
	plt.plot(x,array[:,5],'+',label=file+"_x")
	plt.plot(x,array[:,6],'+',label=file+"_y")

l = plt.legend(prop={'size':6})
l.set_zorder(20)
plt.savefig('gauflux_beam.pdf')
plt.clf()

#files = ['DIP000-190411-2120','DIP000-190223-0052']
files = ['DIP000-190221-0100','DIP000-190223-0030','DIP000-190223-0052','DIP000-190225-0022','DIP000-190225-0044','DIP000-190226-0018','DIP000-190226-0040','DIP000-190227-0036','DIP000-190228-0011','DIP000-190411-2120']
for file in files:
	array = np.loadtxt('output/'+file+'/measurements.txt')
	x = range(0,len(array[array[:,3]==1,6]))
	plt.plot(x,array[array[:,3]==1,6],'+',label=file)

l = plt.legend(prop={'size':6})
l.set_zorder(20)
plt.savefig('skydip.pdf')

plt.clf()
conv_from = [3,4,5,21,22,23,24,25,26]
count = 0
for file in files:
	array = np.loadtxt('output/'+file+'/measurements.txt')
	calval = 0.0
	print(len(array))
	# todelete = np.where(array[:,0] == 23-1)
	# array = np.delete(array, todelete, axis=0)
	# todelete = np.where(array[:,0] == 25-1)
	# array = np.delete(array, todelete, axis=0)
	# todelete = np.where(array[:,0] == 26-1)
	# array = np.delete(array, todelete, axis=0)
	for j in range(1,21):
		for i in range(0,len(array[:,6])):
			try:
				if array[i,3]==j:
					# if array[i,0] == 23 or array[i,0] == 26:
					# 	array = np.delete(array, (i), axis=0)
					# else:
					array[i,0] = conv_from.index(array[i,0])*16 + 4*array[i,1] + array[i,2]
					if array[i,2] == 0:
						calval = array[i,6]
					array[i,6] = array[i,6] / calval
			except:
				continue
		print(len(array))
		count += 1
		plt.plot(array[array[:,3]==j,0],array[array[:,3]==j,6],label=file)

files = ['MOON-190411-1525','CRAB-190410-1529','CRAB-190410-2025','MOON-190223-0528','MOON-190224-0528','MOON-190225-0328']
for file in files:
	array = np.loadtxt('output/'+file+'/_gauflux.txt')
	calval = 0.0
	# print(len(array))
	# todelete = np.where(array[:,0] == 23)
	# array = np.delete(array, todelete, axis=0)
	# todelete = np.where(array[:,0] == 25)
	# array = np.delete(array, todelete, axis=0)
	# todelete = np.where(array[:,0] == 26)
	# array = np.delete(array, todelete, axis=0)
	for i in range(0,len(array[:,4])):
		try:
			# if array[i,0] == 23 or array[i,0] == 26:
			# 	array = np.delete(array, (i), axis=0)
			# else:
			array[i,0] = conv_from.index(array[i,0]-1)*16 + 4*(array[i,1]-1) + (array[i,2]-1)
			if array[i,2] == 1:
				calval = array[i,4]
			array[i,4] = array[i,4] / calval
		except:
			continue
	# print(array[:,4])
	print(len(array))
	# x = range(0,len(array[:,4]))
	plt.plot(array[:,0],array[:,4],label=file)

# l = plt.legend(prop={'size':6})
# l.set_zorder(20)
plt.ylim(0.98,1.02)
plt.savefig('skydip_r.pdf')
plt.clf()
print(count)