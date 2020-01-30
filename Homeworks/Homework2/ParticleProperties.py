# Colin Hauch
# 400B, with Dr. Gurtina Besla
# Homework2


# Location:
# /Users/chauch/Documents/College/Undergraduate/Semester8/400B/400B_hauch/Homeworks/Homework2

# import relavent modules
import numpy as np
import astropy.units as u
from ReadFile import Read


# define a function that takes parameters of information desired and fetches that info from the data file
# used fetched data to return magnitude of the distance in kpc, magnitude of the velocity in km/s, and Mass in units of M⊙
# RETURNS: position, velocity, mass
def ParticleInfo(filename, particleType, particleNumber):

	#First Column type: Particle Type. Type 1 = Dark Matter, Type 2 = Disk Stars, Type 3 = Bulge Stars

	time, total, data = Read(filename)

	index = np.where(data['type'] == particleType)

	# Fetch appropriate variables of a particular particle number within the given particleType and add units
	# NOTE: particleNumber refers to the number of particle OF THE GIVEN particleType
	# 
	# For clarity
	# data 					matrix representation of dataset from filename
	# data['x']				list of all 'x' values from dataset
	# data['x'][index]		list of all 'x' values that have 'type' == particleType
	# data['x'][index][n]	nth element of the list of all 'x' values that have 'type' == particleType
	pos_x = float(data['x'][index][particleNumber])*u.kpc
	pos_y = float(data['y'][index][particleNumber])*u.kpc
	pos_z = float(data['z'][index][particleNumber])*u.kpc

	vel_x = float(data['vx'][index][particleNumber])*(u.km/u.s)
	vel_y = float(data['vy'][index][particleNumber])*(u.km/u.s)
	vel_z = float(data['vz'][index][particleNumber])*(u.km/u.s)

	mass = (float(data['m'][index][particleNumber])* 10**10) *u.Msun


	# caluclate magnitudes (if needed) and round values to three decimal places
	position = np.around(np.sqrt(pos_x**2 + pos_y**2 + pos_z**2), 3)
	velocity = np.around(np.sqrt(vel_x**2 + vel_y**2 + vel_z**2), 3)
	mass = np.around(mass, 3)

	# return values for future use.
	return position, velocity, mass

# test call of the function
Position, Velocity, Mass = ParticleInfo("MW_000.txt",2.0,99)

print(Position, Velocity, Mass)





