# Colin Hauch
# Determine which particles would make good 
# cannidates for view point particle (vpp)

# import revelant modules and classes
import numpy as np
from CenterOfMass import CenterOfMass
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#%matplotlib inline

# my modules
from ReadFile import Read
from MassProfile import MassProfile



# load initial data from M31
# data = np.load("LowResNPY/M31_000.npy")

# declare disk (2) or bulge (1) particles to be used to calc COM
particleType = 2

# find index of particles of only a particlar type
# index = np.where(data['type'] == particleType)

# Determine center of mass of M31
M31com = CenterOfMass("LowRes/M31_000.txt",particleType)
M31_comP = M31com.COM_P(0.1)
M31_comV = M31com.COM_V(M31_comP[0],M31_comP[1],M31_comP[2])

# Determine positions of disk particles relative to COM 
xD = M31com.x - M31_comP[0].value 
yD = M31com.y - M31_comP[1].value 
zD = M31com.z - M31_comP[2].value 

# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)

# Determine velocities of disk particles relatiev to COM motion
vxD = M31com.vx - M31_comV[0].value 
vyD = M31com.vy - M31_comV[1].value 
vzD = M31com.vz - M31_comV[2].value 

# total velocity 
vtot = np.sqrt(vxD**2 + vyD**2 + vzD**2)

# Vectors for r and v 
r = np.array([xD,yD,zD]).T # transposed 
v = np.array([vxD,vyD,vzD]).T


def RotateFrame(posI,velI):
	# input:  3D array of positions and velocities
	# returns: 3D array of rotated positions and velocities such that j is in z direction

	# compute the angular momentum
	L = np.sum(np.cross(posI,velI), axis=0)
	# normalize the vector
	L_norm = L/np.sqrt(np.sum(L**2))


	# Set up rotation matrix to map L_norm to z unit vector (disk in xy-plane)
	
	# z unit vector
	z_norm = np.array([0, 0, 1])
	
	# cross product between L and z
	vv = np.cross(L_norm, z_norm)
	s = np.sqrt(np.sum(vv**2))
	
	# dot product between L and z 
	c = np.dot(L_norm, z_norm)
	
	# rotation matrix
	I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
	v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
	R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

	# Rotate coordinate system
	pos = np.dot(R, posI.T).T
	vel = np.dot(R, velI.T).T
	
	return pos, vel

def rotatedPlot(rnew, vnew):
	# ADD HERE
	# determine the rotated velocity vectors
	# rnew, vnew = RotateFrame(r,v)

	# Rotated M31 Disk - EDGE ON

	# M31 Disk Density 
	fig, ax= plt.subplots(figsize=(10, 10))

	# plot the particle density for M31 , 2D histogram
	#ax.hist2d(rnew[:,0],rnew[:,1], bins=150, norm=LogNorm(), cmap='magma')

	# SCATTER
	ax.scatter(rnew[:,0],rnew[:,1])
	

	#plt.colorbar()

	# Add axis labels
	plt.xlabel('x (kpc)', fontsize=22)
	plt.ylabel('y (kpc)', fontsize=22)

	#set axis limits
	plt.ylim(-50,50)
	plt.xlim(-50,50)

	#adjust tick label font size
	label_size = 22
	matplotlib.rcParams['xtick.labelsize'] = label_size 
	matplotlib.rcParams['ytick.labelsize'] = label_size

	#plt.show()

	return fig, ax


# # determine magnitude of distance of each particle from M31 center of mass
# distance = np.sqrt( (data['x'][index]-comX)**2 + (data['y'][index]-comY)**2 + (data['z'][index]-comZ)**2)



rnew, vnew = RotateFrame(r,v)

fig, ax = rotatedPlot(rnew, vnew)

# np.where((dists >= r) & (dists <= r + dr))
r = 8 # kpc
halfwidth = 0.01 # kpc

distance = np.sqrt( rnew[:,0]**2 + rnew[:,1]**2 + rnew[:,2]**2)

#print(distance)

#disFilter = np.where((distance >= r-halfwidth) & (distance <= r+halfwidth))
disFilter = np.where((distance >= r) & (distance <= r+halfwidth))

print(disFilter[0])


#print(disFilter)
myR = rnew[disFilter]

ax.scatter(myR[:,0],myR[:,1],color="red")


#np.save("cannidates.npy",disFilter[0])


plt.show()












