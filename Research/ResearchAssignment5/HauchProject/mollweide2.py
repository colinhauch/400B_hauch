import matplotlib.pyplot as plt
import numpy as np

# load list of INDEXES of particles that could be used as VPP
cannidateParticles = np.load("cannidates.npy")

# pick one of the indexes
# I picked the 30th particle, as a proof of concept
vppIndex = cannidateParticles[0]

#print(vppIndex)

# load snap data
snap000Data = np.load("LowResNPY/M31_000.npy")
#print(len(snap000Data))

# filter for disk particles (type 2)
index = np.where(snap000Data['type'] == 2)

# fetch the line of data corresponding to VPP
vpp = snap000Data[index][vppIndex]

# print data to make sure I'm not crazy and it works
#print(vpp)

# pull out vpp data to avoid divByZero errors
diskParticleData = np.delete(snap000Data[index],vppIndex)

# fetch all data
# should possibly remove VPP particle...
# but that would change the index of all
# the particles after and numpy seems to 
# handle the divByZero error within the
# trig functions
xdata = diskParticleData['x'] - vpp[2]
ydata = diskParticleData['y'] - vpp[3]
zdata = diskParticleData['z'] - vpp[4]

raBINS = 500
decBINS = 500

# for creating 2D histogram
lon = np.linspace(-np.pi, np.pi,raBINS)
lat = np.linspace(-np.pi/2., np.pi/2.,decBINS)
Lon,Lat = np.meshgrid(lon,lat)

# why convert lat,long to RA,Dec when I can calculate them directly?
RA = np.arctan2(ydata,xdata) #* 180 / np.pi
DEC = np.arcsin(zdata/np.sqrt(xdata**2 + ydata**2 + zdata**2)) #* 180 / np.pi

# print("\nLon (RA)")
# print(lon)
print("\nLat (Dec)")
print(lat)

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

		# count the number of particles in box
		C[i,j] = len(inDECrange[0])

## IT ACTUALLY WORKED I'M GETTING A BEER

		

print(C)



fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')
#arr = np.random.rand(180, 360)
#arr = np.full((180,360))



im = ax.pcolormesh(Lon,Lat,C, cmap=plt.cm.jet)



# print(len(Lon))
# print(len(Lat))

plt.show()












