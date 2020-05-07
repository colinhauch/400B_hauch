# Colin Hauch


# import neccesary modules
import matplotlib.pyplot as plt
import numpy as np
from CenterOfMass2 import CenterOfMass
from matplotlib.colors import LogNorm


def generateImage(snapNumber, vppIndex):
	# PURPOSE: Given a snap number and the address of a particle within the data array
	# generate a single 2D histogram Mollweide plot of the sky from particle specified
	# INPUTS:
	#	snapnumber 	-- string, from "000" to "801"
	#	vppIndex 	-- index of particle designated as the view point particle
	# RETURNS:
	# 	nothing, but saves a plot into a specified folder


	# construct filenames for each galaxy
	filenameM31 = "LowResNPY/M31_"+snapNumber+".npy"
	filenameMW = "LowResNPY/MW_"+snapNumber+".npy"
	filenameM33 = "LowResNPY/M33_"+snapNumber+".npy"


	# load snap data
	snapDataM31 = np.load(filenameM31)
	snapDataMW = np.load(filenameMW)
	snapDataM33 = np.load(filenameM33)

	# collect snap data into a single array
	snapData = np.concatenate((snapDataM31,snapDataMW,snapDataM33))

	# index to filter for disk particles (type 2)
	index = np.where(snapData['type'] == 2)

	# fetch the line of data corresponding to VPP
	vpp = snapData[index][vppIndex]

	# pull out vpp data to avoid divByZero errors
	diskParticleData = np.delete(snapData[index],vppIndex)

	
	### Uncomment following lines to use VPP as view point
	### Otherwise, Center of Mass of M31 will be used
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

	#################################################

	# Declear number of bins for each direction
	raBINS = 60
	decBINS = 60

	# for creating 2D histogram
	# create arrays for each axis, creating the bins
	lon = np.linspace(-np.pi, np.pi,raBINS)
	lat = np.linspace(-np.pi/2., np.pi/2.,decBINS)

	# combine axis arrays to create an array that stores
	# the verticies of each bin
	Lon,Lat = np.meshgrid(lon,lat)

	# Calculate RA and DEC from input data, units of radians
	RA = np.arctan2(ydata,xdata) 
	DEC = np.arcsin(zdata/np.sqrt(xdata**2 + ydata**2 + zdata**2)) 

	# initialize a 2D array to store the values for each bin
	C = np.zeros( ( len(Lon)-1,len(Lat)-1 ) )


	# for each RA slice:
	for i in range(len(lon)-1):

		# determine number of particles within range of RA
		inRArange = np.where((RA>=lon[i])&(RA<lon[i+1]))
		
		# find the declination values of those particles
		remaining = DEC[inRArange]

		# for each declination slice within each RA slice
		for j in range(len(lat)-1):
			
			# find which of those declinations are within a corresponding range
			inDECrange = np.where( (remaining>=lat[j]) & (remaining<lat[j+1]) )

			# count the number of particles in box
			C[i,j] = len(inDECrange[0])

	# Normalize with a logarithmic scale
	# Also transpose, something is backwards with Lat, Lon
	C = np.transpose(np.log(C))

	# Generate figure, as a mollweide plot
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='mollweide')

	# Create title which stores snap number
	# if using VPP, use appropriate title line:
	# plt.title("snap = "+snapNumber+"\nVPP = "+str(vppIndex))
	plt.title("snap = "+snapNumber)

	# Plot the 2D histogram
	im = ax.pcolormesh(Lon,Lat,np.transpose(C), cmap=plt.cm.jet)
	
	### Test plot as a square plot using standard matplotlib functions
	# ax = fig.add_subplot(111)
	# # plt.hist2d(RA,DEC,bins=250, norm=LogNorm(), cmap='magma')
	# plt.hist2d(RA,DEC,bins=250, cmap='magma')

	## Show plot if desired
	#plt.show()

	# save plot in this folder with the snap number string as the name, .png format
	plt.savefig("testSeries1/"+snapNumber)

	# close figure
	plt.close()

	return



# load list of INDEXES of particles that could be used as VPP
cannidateParticles = np.load("cannidates.npy")

# pick one of the indexes
# Selected the 0th particle from list of cannidates
vppIndex = cannidateParticles[0]

# Generate Mollweide Plot for each snap number
for i in range(0,802):

	# gen snap number from 000 to 801
	if len(str(i)) == 3:
		snap = str(i)
	elif len(str(i)) == 2:
		snap = "0" + str(i)
	else:
		snap = "00" + str(i)

	# Call plot generating function
	generateImage(snap,vppIndex)

	# Print snap for to know progression
	print(snap)




