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


			particleMasses = self.m[index][index1]

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

		# From solution:
		# create array of masses within each radius in radii
		Menc = self.MassEnclosed(particleType,radii)

		# calculate circular speed assuming spherical symmetry
		# array of speeds according to "radii", in urnts of kpc/Gyr

        Vcirc = np.round(np.sqrt(G*Menc/radii/u.kpc),2)

        # return array
        return  Vcirc

	def CircularVelocityTotal(self, radii):
		# INPUTS:
		# 		radii -- array with various radius values (kpc)
		#
		# RETURNS:
		# 		array total circular velocities within the given radii
		darkMass, diskMass, bulgeMass = self.MassEnclosedTotal(radii)

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


		return np.around(np.sqrt((G*self.HernquistMass(r,a,darkMass))/r),2)



r = np.arange(0.25, 30.5, 1.15); print(r)



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

plt.xlabel(r'Radius (kpc)', fontsize=22)
plt.ylabel(r'Velocity (km/s)', fontsize=22)

# ax.semilogy(r,M31totalMass,color='red',linewidth=5,label='M31 Total')
# ax.semilogy(r,M31darkMass,color='blue',linewidth=5,label='M31 dark matter')
# ax.semilogy(r,M31diskMass,color='green',linewidth=5,label='M31 disk matter')
# ax.semilogy(r,M31bulgeMass,color='cyan',linewidth=5,label='M31 bulge matter')
# ax.semilogy(r,M31.HernquistMass(r,1,M31darkMass),color='orange',linewidth=5,label='Hernquist Profile, a=1')

# ax.semilogy(r,M31totalMass,color='red',linewidth=5,label='M31 Total')
# ax.semilogy(r,M31darkMass,color='blue',linewidth=5,label='M31 dark matter')
# ax.semilogy(r,M31diskMass,color='green',linewidth=5,label='M31 disk matter')
# ax.semilogy(r,M31bulgeMass,color='cyan',linewidth=5,label='M31 bulge matter')
# ax.semilogy(r,M31.HernquistMass(r,1,M31darkMass),color='orange',linewidth=5,label='Hernquist Profile, a=1')

# ax.semilogy(r,M33totalMass,color='red',linewidth=5,label='M33 Total')
# ax.semilogy(r,M33darkMass,color='blue',linewidth=5,label='M33 dark matter')
# ax.semilogy(r,M33diskMass,color='green',linewidth=5,label='M33 disk matter')
# ax.semilogy(r,M33.HernquistMass(r,1,M33darkMass),color='orange',linewidth=5,label='Hernquist Profile, a=1')


# plt.plot(r,MW.CircularVelocityTotal(r),color='red',linewidth=5,label='MW Total')
# plt.plot(r,MW.HernquistVCirc(r,1,MWdarkMass),color='orange',linewidth=5,label='Hernquist Profile, a=1')

# plt.plot(r,M31.CircularVelocityTotal(r),color='red',linewidth=5,label='M31 Total')
# plt.plot(r,M31.HernquistVCirc(r,1,M31darkMass),color='orange',linewidth=5,label='Hernquist Profile, a=1')

plt.plot(r,M33.CircularVelocityTotal(r),color='red',linewidth=5,label='M33 Total')
plt.plot(r,M33.HernquistVCirc(r,1,M33darkMass),color='orange',linewidth=5,label='Hernquist Profile, a=1')




# Continue plotting for the other redshifts here


# Legend
plt.legend(loc='lower right',fontsize='x-large')
plt.show()









