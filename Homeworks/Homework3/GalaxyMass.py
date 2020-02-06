# Colin Hauch
# 400B, with Dr. Gurtina Besla
# Homework3


# Location:
# /Users/chauch/Documents/College/Undergraduate/Semester8/400B/400B_hauch/Homeworks/Homework3

# import relavent modules
import numpy as np
import astropy.units as u

from ParticleProperties import ParticleInfo
from ReadFile import Read

# Define a function to return the total mass of any desired galaxy component
# particle types are numbers: Halo type (1), Disk type (2), Bulge type (3)
def ComponentMass(filename, particleType):
	# INPUT:	
	#		filename -- the name of the file with the appropriate data
	#		particleType -- Halo type (1), Disk type (2), Bulge type (3)
	#
	# RETURNS:
	#		mass -- the total mass of the specified particle type in units of 10^12 solar masses
	#
	# NOTE: Did not use ParticleProperites.ParticleInfo() and iterate over it because of
	#		variable file length. It seems very difficult to determine the number of particles
	# 		of a particular type within the data file that would be much easier than the
	#		solution I came to below.

	# From ParticleProperties.py
	# Read the data file
	time, total, data = Read(filename)

	# create an a list of index values where the particleType == the type specified
	index = np.where(data['type'] == particleType)

	# Extract a list of the amount of mass over the given index range
	# For clarity
	# data 					matrix representation of dataset from filename
	# data['m']				list of all mass values from dataset
	# data['m'][index]		list of all mass values that have 'type' == particleType
	mass = (data['m'][index]* 10**10) *u.Msun

	# Sum all the entries within the list "mass", then divide 
	# by 10^12 for requested units
	#
	# NOTE: round, then sum them, to help avoid floating point errors
	totalMass = sum(mass)/(10**12)

	# round value to requested precision (3 points after decimal)
	totalMass = np.around(totalMass, 3)

	# return the total mass of the specified particle type
	return totalMass


# Printing for homwork
# 
# print(ComponentMass("MW_000.txt",1))
# print(ComponentMass("MW_000.txt",2))
# print(ComponentMass("MW_000.txt",3))
# print("\n")
# print(ComponentMass("M31_000.txt",1))
# print(ComponentMass("M31_000.txt",2))
# print(ComponentMass("M31_000.txt",3))
# print("\n")
# print(ComponentMass("M33_000.txt",1))
# print(ComponentMass("M33_000.txt",2))








