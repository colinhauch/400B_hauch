# Colin Hauch




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

class GalaxyOrbit:
	# create a class that is utilized to calculate the orbit of 
	# any of the three galaxies.

	def __init__(self, filenameFormat):
		# INPUTS:
		#
		# 	filenameFormat -- string, format of the
		#					  name of the file to store integrated orbit
		
		### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### store the output file name
        self.filename = filenameFormat
        

        # get the current position and velocity of each galaxy
        MW_centerOfMass  = CenterOfMass("MW_000.txt",2)
        M31_centerOfMass = CenterOfMass("M31_000.txt",2)
        M33_centerOfMass = CenterOfMass("M33_000.txt",2)

        # get the position center of mass for each galaxy
        MW_comP  = MW_centerOfMass.COM_P(0.1)
        M31_comP = M31_centerOfMass.COM_P(0.1)
        M33_comP = M33_centerOfMass.COM_P(0.1)

        # get the velocity center of mass for each galaxy
        # TODO: check if phrase "velocity center of mass" actually means something
        MW_comV  = MW_centerOfMass.COM_V(MW_comP[0], MW_comP[1], MW_comP[2])
        M31_comV = M31_centerOfMass.COM_V(M31_comP[0], M31_comP[1], M31_comP[2])
        M33_comV = M33_centerOfMass.COM_V(M33_comP[0], M33_comP[1], M33_comP[2])

        # store difference vectors between each galaxy vector type (pos/vel)
        # may have to be more careful about relative to which galaxy these 
        # calculations are done. In the mean time, all six vectors are determined
        # key: MW_M31 means "MW - M31"

        # Position
        self.r0__MW_M31  = (MW_comP - M31_comP).value
        self.r0__MW_M33  = (MW_comP - M33_comP).value

        self.r0__M31_MW  = (M31_comP - MW_comP).value
        self.r0__M31_M33 = (M31_comP - M33_comP).value

        self.r0__M33_M31 = (M33_comP - M31_comP).value
        self.r0__M33_MW  = (M33_comP - MW_comP).value

        # Velocity
        self.v0__MW_M31  = (MW_comV - M31_comV).value
        self.v0__MW_M33  = (MW_comV - M33_comV).value

        self.v0__M31_MW  = (M31_comV - MW_comV).value
        self.v0__M31_M33 = (M31_comV - M33_comV).value

        self.v0__M33_M31 = (M33_comV - M31_comV).value
        self.v0__M33_MW  = (M33_comV - MW_comV).value

        # Get scale lengths for each mass component (using values for HW7, M31)
        self.MW_rdisk 	= 5
        self.M31_rdisk 	= 5
        self.M33_rdisk 	= 5

        self.MW_rbulge 	= 5
        self.M31_rbulge = 5
        self.M33_rbulge = 5

        self.MW_rhalo 	= 5
        self.M31_rhalo 	= 5
        self.M33_rhalo 	= 5

        # Get the mass of each component in each galaxy
        # in units of M_sun
        self.MW_Mhalo 	= ComponentMass("MW_000.txt",1)*1e12
        self.M31_Mhalo 	= ComponentMass("M31_000.txt",1)*1e12
        self.M33_Mhalo 	= ComponentMass("M33_000.txt",1)*1e12

        self.MW_Mdisk 	= ComponentMass("MW_000.txt",2)*1e12
        self.M31_Mdisk 	= ComponentMass("M31_000.txt",2)*1e12
        self.M33_Mdisk 	= ComponentMass("M33_000.txt",2)*1e12

        self.MW_Mbulge 	= ComponentMass("MW_000.txt",3)*1e12
        self.M31_Mbulge = ComponentMass("M31_000.txt",3)*1e12
        self.M33_Mbulge = ComponentMass("M33_000.txt",3)*1e12

    def comp_mag(vector):
    	# INPUT:
    	#
    	# 	vector -- an array with exactly three components
    	#
    	# RETURNS:
    	#
    	#	magnitude of vector
        return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2) 

	def HernquistAccel(self, M, r_a, posVEC): 
        # it is easiest if you take as an input the position VECTOR 
        # Function that determines the acceleration of a hernquist profile
        # INPUTS:
        #   M -- total mass of halo or bulge
        #   r_a -- scale lengnth corresponding to halo or bulge
        #   posVEC -- position vector 
        #
        # RETURNS:
        #   acceleration vector for hernquist profile
        
        # Store the magnitude of the position vector
        rmag = comp_mag(posVEC)
        
        # Store the Acceleration
        Hern = -((self.G*M) / (rmag*(r_a + rmag)**2)) * posVEC

        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern

	def MiyamotoNagaiAccel(self, M, rd, posVEC):
        # it is easiest if you take as an input a position VECTOR  r 
        # Function that determines the acceleration of a disk particles
        # using a profile better suited for them
        # INPUTS:
        #   M -- total mass of halo or bulge
        #   rd -- disk radius of M31
        #   posVEC -- position vector 
        #
        # RETURNS:
        #   acceleration vector

        zd = rd / 5.0
        R = np.sqrt(posVEC[0]**2 + posVEC[1]**2)
        B = rd + np.sqrt(posVEC[2]**2 + zd**2)

        ZSTUFF = B/np.sqrt(posVEC[2]**2 + zd**2)
        zModification = np.array([1,1,ZSTUFF]) 

        # the zModification allows for a different value for the z component of the acceleration
        return ((-self.G * M)/(R**2 + B**2)**1.5)*posVEC*zModification

    def Galaxy_Accel(self, charaVEC, posVEC):
    	# Function that takes a position vector and sums all 
        # acceleration vectors from each galaxy component
        #
        # INPUTS:
        #	charaVEC -- short for "characteristic vector": vector containing
        #					the mass values of each component of the galaxy and
        # 					the scale lengths of each component.
        #					[Mhalo,Mdisk,Mbulge,rhalo,rdisk,rbulge]
        #   posVEC -- position vector 
        #
        # RETURNS:
        #   summed acceleration vector

        # "unpack" input characteristic vector
        Mhalo  = charaVEC[0]
        Mdisk  = charaVEC[1]
        Mbulge = charaVEC[2]

        rhalo  = charaVEC[3]
        rdisk  = charaVEC[4]
        rbulge = charaVEC[5]


        # Call the previous functions for the halo, bulge and disk
        # these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        haloAcceleration = self.HernquistAccel(Mhalo, rhalo, posVEC)
        diskAcceleration = self.MiyamotoNagaiAccel(Mdisk, rdisk, posVEC)
        bulgeAcceleration = self.HernquistAccel(Mbulge, rbulge, posVEC)
            
        # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return np.sum([haloAcceleration,diskAcceleration,bulgeAcceleration], axis=0)
    
	def LeapFrog(self, galaxyName, dt, posVEC, velVEC): 
        # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        # Integrates to solve for the orbit of ONE galaxy
        # INPUTS:
        #	
        #	galaxyName -- string of galaxy name ("MW","M31","M33")
        #   dt -- time interval for integration
        #   posVEC -- position vector 
        #   velVEC -- velocity vector
        #
        # RETURNS:
        #   new position and velocity vectors, packed

        # characteristicVector -- short for "characteristic vector": 
        # vector containing the mass values of each component of the 
        # galaxy and the scale lengths of each component. 
        # [Mhalo,Mdisk,Mbulge,rhalo,rdisk,rbulge]
        characteristicVector = np.zeros(6)

        if galaxyName == "MW":
        	characteristicVector[0] = self.MW_Mhalo
        	characteristicVector[1] = self.MW_Mdisk
        	characteristicVector[2] = self.MW_Mbulge

        	characteristicVector[3] = self.MW_rhalo
        	characteristicVector[4] = self.MW_rdisk
        	characteristicVector[5] = self.MW_rbulge

        else if galaxyName == "M31":
        	characteristicVector[0] = self.M31_Mhalo
        	characteristicVector[1] = self.M31_Mdisk
        	characteristicVector[2] = self.M31_Mbulge

        	characteristicVector[3] = self.M31_rhalo
        	characteristicVector[4] = self.M31_rdisk
        	characteristicVector[5] = self.M31_rbulge

    	else if galaxyName == "M33":
        	characteristicVector[0] = self.M33_Mhalo
        	characteristicVector[1] = self.M33_Mdisk
        	characteristicVector[2] = self.M33_Mbulge

        	characteristicVector[3] = self.M33_rhalo
        	characteristicVector[4] = self.M33_rdisk
        	characteristicVector[5] = self.M33_rbulge
        else:
        	print("ERROR: galaxy name not specified in func: LeapFrog()")

        
        # predict the position at the next half timestep
        rhalf = posVEC + (velVEC*dt/2)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        accel = self.Galaxy_Accel(characteristicVector,rhalf) # acceleration of M31 at the half step
        vnew = velVEC + (accel*dt)
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don"t know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = posVEC + dt*(velVEC + vnew)/2
        
        # **** return the new position and velcoity vectors
        return rnew,vnew

    def OrbitIntegration(self, t0, dt, tmax):
    	# Loops over LeapFrog integrator to solve the equations 
    	# of motion and compute the future orbit of all three
    	# galaxies for tmax Gyr into the future
        # INPUTS:
        #   t0 -- starting time
        #   dt -- time interval 
        #   tmax -- final time
        #
        # RETURNS:
        #   nothing?

        # initialize the time to the input starting time
        t = t0

		# initialize the postion and velocity vectors to their starting values
        MW_pos  = self.r0__MW_M31
        MW_vel  = self.v0__MW_M33

        M31_pos = self.r0__M31_MW
        M31_vel = self.v0__M31_M33

        M33_pos = self.r0__M33_MW
        M33_vel = self.v0__M33_M31

        # initialize empty arrays of size :  
        # (rows int(tmax/dt)+1 , columns 7)
        MW_orbit  = np.zeros( (np.int(tmax/dt)+1, 7) )
        M31_orbit = np.zeros( (np.int(tmax/dt)+1, 7) )
        M33_orbit = np.zeros( (np.int(tmax/dt)+1, 7) )

        # initialize the first row of the orbits
        MW_orbit[0]  = t0, *tuple(MW_pos), *tuple(MW_vel)
        M31_orbit[0] = t0, *tuple(M31_pos), *tuple(M31_vel)
        M33_orbit[0] = t0, *tuple(M33_pos), *tuple(M33_vel)

        # initialize a counter for the orbit. 
        # since we already set the 0th values, we start the counter at 1
        i = 1 

        # start the integration (advancing in time steps and 
        # computing LeapFrog at each step)
        #
        # as long as t has not exceeded the maximal time:
        while ( t < tmax ):  

        	# advance the time by one timestep, dt
            t += dt

            # advance the position and velocity of each galaxy
            # using the LeapFrog scheme
            MW_pos, MW_vel   = self.LeapFrog("MW", dt, MW_pos, MW_vel)
            M31_pos, M31_vel = self.LeapFrog("MW", dt, M31_pos, M31_vel)
            M33_pos, M33_vel = self.LeapFrog("MW", dt, M33_pos, M33_vel)

            # store the new states in the first column of the ith row
            # of the appropriate vector
            MW_orbit[i]  = t, *tuple(MW_pos), *tuple(MW_vel)
            M31_orbit[i] = t, *tuple(M31_pos), *tuple(M31_vel)
            M33_orbit[i] = t, *tuple(M33_pos), *tuple(M33_vel)

            # update counter i, where i is keeping track of 
            # the number of rows (i.e. the number of time steps)
            i += 1

        MW_filename = "MW_" + self.filename
        M31_filename = "M31_" + self.filename
        M33_filename = "M33_" + self.filename


        # write the data to the appropriate file
        # MW file
        np.savetxt(MW_filename, MW_orbit, fmt = "%11.3f"*7, comments="#", 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format("t", "x", "y", "z", "vx", "vy", "vz"))
        
        # M31 file
        np.savetxt(M31_filename, M31_orbit, fmt = "%11.3f"*7, comments="#", 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format("t", "x", "y", "z", "vx", "vy", "vz"))
        
        # M33 file
        np.savetxt(M33_filename, M33_orbit, fmt = "%11.3f"*7, comments="#", 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format("t", "x", "y", "z", "vx", "vy", "vz"))




