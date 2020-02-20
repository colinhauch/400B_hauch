# Colin Hauch
# 400B, with Dr. Gurtina Besla
# Homework5
# Conferred with Sam Andrews, Collin Lewin, and Chirag Rathi


# Location:
# /Users/chauch/Documents/College/Undergraduate/Semester8/400B/400B_hauch/Homeworks/Homework5

# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass

G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

# a class that will help us analyze the mass profile given a galaxy and a snap number
class MassProfile:
	
	# define the initialization function for the class
	def __init__(self, galaxy, snap):
		# INPUTS:
		# 		galaxy -- a string with Galazy Name ("MW", "M31", "M33")
		# 		snap -- snapshot number as an integer

		# Create a string to save as the filename
		# add a string of the filenumber to the value “000” 
		ilbl = "000" + str(snap)
		# remove all but the last 3 digits
		ilbl = ilbl[-3:]
		self.filename="%s_"%(galaxy) + ilbl + ".txt"

		# read data in the given file using Read
		self.time, self.total, self.data = Read(self.filename)

		# store the mass, positions, and velocities
		self.m = self.data['m']
		self.x = self.data['x']
		self.y = self.data['y']
		self.z = self.data['z']
		self.vx = self.data['vx']
		self.vy = self.data['vy']
		self.vz = self.data['vz']

		# store the string of the name of the galaxy of interest
		self.gname = galaxy

	def MassEnclosed(self, particleType, radii):
		# INPUTS:
		#		particleType -- Particle Type. Type 1 = Dark Matter, Type 2 = Disk Stars, Type 3 = Bulge Stars
		# 		radii -- array with various radius values (kpc)
		#
		# RETURNS:
		# 		array of masses (Msun)

		# create center of mass object
		COM = CenterOfMass(self.filename, particleType)

		# establish coordinates of the center of mass positon
		COM_P = COM.COM_P(0.1)

		# calculate magnitude of COM_P 
		comRadius = np.around(np.linalg.norm(COM_P),3)

		index = np.where(self.data['type'] == particleType)

		xdists = self.x[index] - COM_P[0].value
		ydists = self.y[index] - COM_P[1].value
		zdists = self.z[index] - COM_P[2].value

		particleRadii = np.sqrt(xdists**2 + ydists**2 + zdists**2)

		# create an empty list 
		massValues = np.array([])

		for element in radii:
			# for each distance given, find all the particles within that distance
			index1 = np.where(particleRadii <= element)

			particleMasses = self.m[index1]

			massWithinElement = np.sum(particleMasses)
			#print(massWithinElement)

			massValues = np.append(massValues, massWithinElement)

		return massValues*1e10*u.Msun
		
	def MassEnclosedTotal(self, radii):

		if self.gname == "M33":
			darkMass  = self.MassEnclosed(1,radii)
			diskMass  = self.MassEnclosed(2,radii)
			bulgeMass = 0 * u.Msun
		else:
			darkMass  = self.MassEnclosed(1,radii)
			diskMass  = self.MassEnclosed(2,radii)
			bulgeMass = self.MassEnclosed(3,radii)

		return darkMass + diskMass + bulgeMass

	def HernquistMass(self,r,a,darkMass):

		return darkMass * r**2 / (a+r)**2

	def CircularVelocity(self, particleType, radii):

		index = np.where(self.data['type'] == particleType)



		xdists = self.x[index] - COM_P[0].value
		ydists = self.y[index] - COM_P[1].value
		zdists = self.z[index] - COM_P[2].value

		particleRadii = np.sqrt(xdists**2 + ydists**2 + zdists**2)

		circVelos = np.array([])

		for element in particleRadii:

			M = MassEnclosed([element])

			circVelo = np.sqrt(G*M/element)

			circVelos = np.append(circVelos, circVelo)

		return circVelos



r = np.arange(0.25, 30.5, 1.5); print(r)

MW = MassProfile("MW", 0)
print(MW.MassEnclosed(1,r))
#print(test1.MassEnclosed(1,2))







