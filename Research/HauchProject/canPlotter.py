# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
#%matplotlib inline

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile



# Code for plotting contours
#from https://gist.github.com/adrn/3993992

#import scipy.optimize as so

# Create a COM of object for M31 Disk Using Code from Assignment 4
COMD = CenterOfMass("LowRes/M31_000.txt",2)

# Compute COM of M31 using disk particles
COMP = COMD.COM_P(0.1)
COMV = COMD.COM_V(COMP[0],COMP[1],COMP[2])

# Determine positions of disk particles relative to COM 
xD = COMD.x - COMP[0].value 
yD = COMD.y - COMP[1].value 
zD = COMD.z - COMP[2].value 

# total magnitude
rtot = np.sqrt(xD**2 + yD**2 + zD**2)

# Determine velocities of disk particles relatiev to COM motion
vxD = COMD.vx - COMV[0].value 
vyD = COMD.vy - COMV[1].value 
vzD = COMD.vz - COMV[2].value 

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

def rotatedDensityPlot():
	# ADD HERE
	# determine the rotated velocity vectors
	rnew, vnew = RotateFrame(r,v)

	# Rotated M31 Disk - EDGE ON

	# M31 Disk Density 
	fig, ax= plt.subplots(figsize=(10, 10))

	# plot the particle density for M31 , 2D histogram
	# ADD HERE
	ax.hist2d(rnew[:,0],rnew[:,1], bins=150, norm=LogNorm(), cmap='magma')

	#plt.colorbar()

	# Add axis labels
	plt.xlabel('x (kpc)', fontsize=22)
	plt.ylabel('y (kpc)', fontsize=22)

	#set axis limits
	plt.ylim(-40,40)
	plt.xlim(-40,40)

	#adjust tick label font size
	label_size = 22
	matplotlib.rcParams['xtick.labelsize'] = label_size 
	matplotlib.rcParams['ytick.labelsize'] = label_size

	#plt.show()

	return fig, ax


#fig, ax = rotatedDensityPlot()














#densityPlot()
