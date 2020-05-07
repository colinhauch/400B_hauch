import matplotlib.pyplot as plt
import numpy as np
from CenterOfMass import CenterOfMass
from matplotlib.colors import LogNorm

def generateImage(snapNumber, vppIndex):

	# snapNumber is a string!!

	filenameM31 = "LowResNPY/M31_"+snapNumber+".npy"
	filenameMW = "LowResNPY/MW_"+snapNumber+".npy"
	filenameM33 = "LowResNPY/M33_"+snapNumber+".npy"


	# load snap data
	snapDataM31 = np.load(filenameM31)
	snapDataMW = np.load(filenameMW)
	snapDataM33 = np.load(filenameM33)
	#print(len(snap000Data))

	snapData = np.concatenate((snapDataM31,snapDataMW,snapDataM33))
	

	# filter for disk particles (type 2)
	index = np.where(snapData['type'] == 2)

	# fetch the line of data corresponding to VPP
	vpp = snapData[index][vppIndex]



	# print data to make sure I'm not crazy and it works
	#print(vpp)

	# pull out vpp data to avoid divByZero errors
	diskParticleData = np.delete(snapData[index],vppIndex)

	#print(len(diskParticleData))

	# fetch all data
	# xdata = diskParticleData['x'] - vpp[2]
	# ydata = diskParticleData['y'] - vpp[3]
	# zdata = diskParticleData['z'] - vpp[4]

	############## Center of Mass of M31 ############

	M31com = CenterOfMass(filenameM31,2)
	M31_comP = M31com.COM_P(0.1)

	xdata = diskParticleData['x'] - M31_comP[0].value 
	ydata = diskParticleData['y'] - M31_comP[1].value 
	zdata = diskParticleData['z'] - M31_comP[2].value 

	##########################

	raBINS = 60
	decBINS = 60

	# for creating 2D histogram
	lon = np.linspace(-np.pi, np.pi,raBINS)
	lat = np.linspace(-np.pi/2., np.pi/2.,decBINS)
	Lon,Lat = np.meshgrid(lon,lat)

	# why convert lat,long to RA,Dec when I can calculate them directly?
	RA = np.arctan2(ydata,xdata) #* 180 / np.pi
	DEC = np.arcsin(zdata/np.sqrt(xdata**2 + ydata**2 + zdata**2)) #* 180 / np.pi

	# print("\nLon (RA)")
	# print(lon)
	# print("\nLat (Dec)")
	# print(lat)

	C = np.zeros( ( len(Lon)-1,len(Lat)-1 ) )

	#print(len(RA))

	for i in range(len(lon)-1):

		# determine number of particles within range of RA
		inRArange = np.where((RA>=lon[i])&(RA<lon[i+1]))
		
		#print("length of inRArange = "+str(len(inRArange[0])))

		# find the declination values of those particles
		remaining = DEC[inRArange]

		for j in range(len(lat)-1):
			
			# find which of those declinations are within a corresponding range
			inDECrange = np.where( (remaining>=lat[j]) & (remaining<lat[j+1]) )

			#print(len(inDECrange[0]))

			# count the number of particles in box
			C[i,j] = len(inDECrange[0])

	## IT ACTUALLY WORKED I'M GETTING A BEER

	C = np.log(C)
	# print(C)
	# print(np.amax(C))

	# print(C)



	fig = plt.figure()
	ax = fig.add_subplot(111, projection='mollweide')
	plt.title("snap = "+snapNumber+"\nVPP = "+str(vppIndex))
	im = ax.pcolormesh(Lon,Lat,C, cmap=plt.cm.jet)
	
	
	# ax = fig.add_subplot(111)
	# # plt.hist2d(RA,DEC,bins=250, norm=LogNorm(), cmap='magma')
	# plt.hist2d(RA,DEC,bins=250, cmap='magma')


	#plt.show()

	plt.savefig("LowResSeries3/"+snapNumber)

	plt.close()

	return



# load list of INDEXES of particles that could be used as VPP
cannidateParticles = np.load("cannidates.npy")

# pick one of the indexes
# I picked the 30th particle, as a proof of concept
vppIndex = cannidateParticles[0]




for i in range(0,802):

	# gen snap number from 000 to 801
	if len(str(i)) == 3:
		snap = str(i)
	elif len(str(i)) == 2:
		snap = "0" + str(i)
	else:
		snap = "00" + str(i)


	generateImage(snap,vppIndex)
	print(snap)




