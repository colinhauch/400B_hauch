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
import matplotlib.pyplot as plt
import matplotlib

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
		# INPUTS:
		# 		radii -- array with various radius values (kpc)
		#
		# RETURNS:
		# 		array of total mass enclosed within the radii (Msun)

		# if the galaxy is M33 (which doesn't have a bulge), make an exception
		if self.gname == "M33":

			# determine mass of each component but not the bulge
			darkMass  = self.MassEnclosed(1,radii)
			diskMass  = self.MassEnclosed(2,radii)
			bulgeMass = 0 * u.Msun

		# otherwise,
		else:
			# determine mass of each component
			darkMass  = self.MassEnclosed(1,radii)
			diskMass  = self.MassEnclosed(2,radii)
			bulgeMass = self.MassEnclosed(3,radii)

		# return the sum of the masses
		return darkMass, diskMass, bulgeMass

	def HernquistMass(self,r,a,darkMass):
		# INPUTS:
		#		r -- radius
		# 		a -- scale length used in this profile
		#		darkMass -- amount of dark matter to be used in the profile
		#
		# RETURNS:
		# 		computed mass enclosed according to the profile

		return darkMass * r**2 / (a+r)**2

	def CircularVelocity(self, particleType, radii):
		# INPUTS:
		#		particleType -- Particle Type. Type 1 = Dark Matter, Type 2 = Disk Stars, Type 3 = Bulge Stars
		# 		radii -- array with various radius values (kpc)
		#
		# RETURNS:
		# 		array circular velocities

		# create an index of the particles of the desired type
		index = np.where(self.data['type'] == particleType)

		# Subtract the COM position from each component value for every particle
		# of the desired type
		xdists = self.x[index] - COM_P[0].value
		ydists = self.y[index] - COM_P[1].value
		zdists = self.z[index] - COM_P[2].value

		# calculate magnitude of position vectors
		particleRadii = np.sqrt(xdists**2 + ydists**2 + zdists**2)

		# create an empty array
		circVelos = np.array([])

		# for each radii in particle radii list, 
		for element in particleRadii:

			# pass that single value radii, as an array, to the MassEnclosed function
			# which returns an ARRAY. The returned array in this case will be of length 1
			M = MassEnclosed([element])

			# calculate the circular velocity for that radii
			circVelo = np.sqrt(G*M[0]/element)

			# append the calculated circular velocity to the array that was initlized before
			circVelos = np.append(circVelos, circVelo)

		return circVelos

	def CircularVelocityTotal(self, radii):
		# INPUTS:
		# 		radii -- array with various radius values (kpc)
		#
		# RETURNS:
		# 		array total circular velocities within the given radii
		darkMass, diskMass, bulgeMass = MassEnclosedTotal(radii)

		totalMass = darkMass + diskMass + bulgeMass

		return np.sqrt((G*totalMass)/radii)

	def HernquistVCirc(self,r,a,darkMass):
		# INPUTS:
		#		radii -- 
		# 		a -- scale length used in this profile
		#		darkMass -- amount of dark matter to be used in the profile
		#
		# RETURNS:
		# 		computed circular velocity according to the hernquist profile mass
		# 		within the given radius


		return np.around(np.sqrt((G*HernquistMass(radii,a,darkMass))/radii),2)



r = np.arange(0.25, 30.5, 1.5); print(r)



MW = MassProfile("MW", 0)
M31 = MassProfile("M31", 0)
M33 = MassProfile("M33", 0)

MWdarkMass, MWdiskMass, MWbulgeMass = MW.MassEnclosedTotal(r)
M31darkMass, M31diskMass, M31bulgeMass = M31.MassEnclosedTotal(r)
M33darkMass, M33diskMass, M33bulgeMass = M33.MassEnclosedTotal(r)

MWtotalMass = MWdarkMass + MWdiskMass + MWbulgeMass
M31totalMass = M31darkMass + M31diskMass + M31bulgeMass
M33totalMass = M33darkMass + M33diskMass + M33bulgeMass


fig,ax = plt.subplots(figsize=(10,8))


#adjust tick label font size
label_size = 22

matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# 
plt.plot(r, MWdarkMass,linewidth = 5, label='MW dark matter',color="blue",linestyle="-")
plt.plot(r, MWdiskMass,linewidth = 5, label='MW disk matter',color="blue",linestyle=":")
plt.plot(r, MWbulgeMass,linewidth = 5, label='MW bulge matter',color="blue",linestyle="--")

# 
plt.plot(r, M31darkMass,linewidth = 5, label='M31 dark matter',color="red",linestyle="-")
plt.plot(r, M31diskMass,linewidth = 5, label='M31 disk matter',color="red",linestyle=":")
plt.plot(r, M31bulgeMass,linewidth = 5, label='M31 bulge matter',color="red",linestyle="--")

# 
plt.plot(r, M33darkMass,linewidth = 5, label='M33 dark matter',color="green",linestyle="-")
plt.plot(r, M33diskMass,linewidth = 5, label='M33 disk matter',color="green",linestyle=":")





# Continue plotting for the other redshifts here


# Legend
plt.legend(loc='lower right',fontsize='x-large')
plt.show()









