# Colin Hauch



# import relavent modules
import numpy as np
import astropy.units as u
from ReadFile import Read

def findParticlesOfInterest(filename,particleType=2):
	# FUNCTION:
	# INPUTS:
	#	filename -- appropriate file for MW, M31, or M33
	# 	particleType -- 1 (halo), 2 (disk), 3 (bulge)
	#					defaults to disk particles, because the sun
	#					is considered a disk particle
	#
	# RETURNS:
	#	list of the address of the particles of interest
	# 	within galaxy

	time, total, data = Read(filename)
	index = np.where(data['type'] == particleType)


	# "POI" -- particles of interest
	POI = np.array([])

	# percent range of similarity to sun
	# EX: "massSIM = 0.1" means within 10% of the sun's mass
	posSIM = 0.1
	velSIM = 0.1
	massSIM = 0.1

	#sunMass 	= u.Msun.value
	sunRadius 	= 8 * u.kpc
	#sunVelocity = 0 * u.km / u.s

	print(len(index))

	for i in range(len(index)):
		# Fetch appropriate variables of a particular particle number within the given particleType and add units
		# NOTE: i refers to the number of particle OF THE GIVEN particleType
		# 
		# For clarity
		# data 					matrix representation of dataset from filename
		# data['x']				list of all 'x' values from dataset
		# data['x'][index]		list of all 'x' values that have 'type' == particleType
		# data['x'][index][n]	nth element of the list of all 'x' values that have 'type' == particleType
		pos_x = float(data['x'][index][i])*u.kpc
		pos_y = float(data['y'][index][i])*u.kpc
		pos_z = float(data['z'][index][i])*u.kpc

		vel_x = float(data['vx'][index][i])*(u.km/u.s)
		vel_y = float(data['vy'][index][i])*(u.km/u.s)
		vel_z = float(data['vz'][index][i])*(u.km/u.s)

		mass = (float(data['m'][index][i])* 10**10) *u.Msun


		# caluclate magnitudes (if needed) and round values to three decimal places
		position = np.around(np.sqrt(pos_x**2 + pos_y**2 + pos_z**2), 3)
		velocity = np.around(np.sqrt(vel_x**2 + vel_y**2 + vel_z**2), 3)
		mass = np.around(mass, 3)

		if (sunRadius*(1-posSIM)) < abs(position) < (sunRadius*(1+posSIM)):
			#if (sunMass*(1-massSIM)) < mass < (sunMass*(1+massSIM)):
			#if velocity
			POI = np.append(POI,i)
			print("found one: particle #"+str(i))
		else:
			break #print("nope")

	return POI


			



	# test call of the function
	# Position, Velocity, Mass = ParticleInfo("MW_000.txt",2.0,99)

	# print(Position, Velocity, Mass)



MW_POI = findParticlesOfInterest("MW_000.txt")

print(len(MW_POI))
print(MW_POI)