# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# look for "****"

# Colin Hauch



# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
#from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass


# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass


class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self,filename): 
        # INPUTS:
        #
        # filename -- string, name of the file to store integrated orbit

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33
        M33_com = CenterOfMass("M33_000.txt",2)
        # Store postion and velocity vectors
        M33_comP = M33_com.COM_P(0.1)
        M33_comV = M33_com.COM_V(M33_comP[0], M33_comP[1], M33_comP[2])
        
   
        
        ### get the current pos/vel of M31 
        M31_com = CenterOfMass("M31_000.txt",2)
        M31_comP = M31_com.COM_P(0.1)
        M31_comV = M31_com.COM_V(M31_comP[0], M31_comP[1], M31_comP[2])


        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        self.r0 = (M33_comP - M31_comP).value
        self.v0 = (M33_comV - M31_comV).value
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5 

        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass("M31_000.txt",2)*1e12 #in M_sun
        
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1

        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass("M31_000.txt",3)*1e12
        
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62

        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass("M31_000.txt",1)*1e12
     
    
    
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
        
        ### **** Store the magnitude of the position vector
        rmag = comp_mag(posVEC)
        
        ### *** Store the Acceleration
        Hern = -((self.G*M) / (rmag*(r_a + rmag)**2)) * posVEC

        #follow the formula in the HW instructions
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

        zd = self.rdisk / 5.0
        R = np.sqrt(posVEC[0]**2 + posVEC[1]**2)
        B = rd + np.sqrt(posVEC[2]**2 + zd**2)

        ZSTUFF = B/np.sqrt(posVEC[2]**2 + zd**2)
        zModification = np.array([1,1,ZSTUFF]) 

        # the zModification allows for a different value for the z component of the acceleration
        return ((-self.G * M)/(R**2 + B**2)**1.5)*posVEC*zModification

        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        
     
    
    def M31Accel(self, posVEC): 
        # Function that takes a position vector and sums all 
        # acceleration vectors from each galaxy component
        #
        # INPUTS:
        #   posVEC -- position vector 
        #
        # RETURNS:
        #   summed acceleration vector
        

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        haloAcceleration = self.HernquistAccel(self.Mhalo, self.rhalo, posVEC)
        diskAcceleration = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, posVEC)
        bulgeAcceleration = self.HernquistAccel(self.Mbulge, self.rbulge, posVEC)
            
        # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return np.sum([haloAcceleration,diskAcceleration,bulgeAcceleration], axis=0)
    
    
    
    def LeapFrog(self, dt, posVEC, velVEC): 
        # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        # Integrates to solve for the orbit of M33
        # INPUTS:
        #   dt -- time interval for integration
        #   posVEC -- position vector 
        #   velVEC -- velocity vector
        #
        # RETURNS:
        #   new position and velocity vectors, packed
        
        # predict the position at the next half timestep
        rhalf = posVEC + (velVEC*dt/2)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        accel = self.M31Accel(rhalf) # acceleration of M31 at the half step
        vnew = velVEC + (accel*dt)
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don"t know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = posVEC + dt*(velVEC + vnew)/2
        
        # **** return the new position and velcoity vectors
        return rnew,vnew
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        # Loops over LeapFrog integrator to solve the equations of motion and compute the 
        # future orbit of M33 for 10 Gyr into the future
        # INPUTS:
        #   t0 -- starting time
        #   dt -- time interval 
        #   tmax -- final time
        #
        # RETURNS:
        #   

        # initialize the time to the input starting time
        t = t0

        # initialize the postion and velocity vectors to their starting values
        pos = self.r0
        vel = self.v0
        
        # initialize an empty array of size :  rows int(tmax/dt)+1  , columns 7
        orbit = np.zeros( (np.int(tmax/dt)+1, 7) )
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # orbit[0] = t0, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while ( t < tmax ):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t += dt
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            pos, vel = self.LeapFrog(dt, pos, vel)

            # **** store the new time in the first column of the ith row
            orbit[i] = t, *tuple(pos), *tuple(vel)
            #orbit[i] = t, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
         
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments="#", 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format("t", "x", "y", "z", "vx", "vy", "vz"))
        
        # there is no return function
    def comp_mag(vector):
        return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2) 

