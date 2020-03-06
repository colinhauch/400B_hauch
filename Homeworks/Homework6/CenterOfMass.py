# Homework 4
# Center of Mass Position and Velocity
# Colin Hauch

###############################################################################
# Keep in mind this is just a template; you don't need to follow every step and feel free to change anything.
# We also strongly encourage you to try to develop your own method to solve the homework.
###############################################################################
# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read


class CenterOfMass:
    # Class to define COM position and velocity properties of a given galaxy 
    # and simulation snapshot
    
    
    def __init__(self, filename, ptype):
        # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]

        # write your own code to complete this for positions and velocities
        self.x  = self.data['x'][self.index]
        self.y  = self.data['y'][self.index]
        self.z  = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
        # Function to compute the center of mass position or velocity generically
        # input: array (a,b,c) of positions or velocities and the mass
        # returns: 3 floats  (the center of mass coordinates)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = np.dot(a,m)/np.sum(m)
        # ycomponent Center of mass
        Bcom = np.dot(b,m)/np.sum(m)
        # zcomponent Center of mass
        Ccom = np.dot(c,m)/np.sum(m)
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta, VolDec):
        # Function to specifically return the center of mass position and velocity                                         
        # input:                                                                                                                                                                                        
        #       delta (tolerance)   
        #       VolDec - defines the amount by which RMAX is decreased                                                                                      
        # returns: One vector, with rows indicating:                                                                                                                                                                            
        #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)  
        # MODIFICATION: division by 2 changed to division by VolDec                                                         
        RMAX = max(RNEW)/VolDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW <= RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/VolDec
            # check this.                                                                                              
            #print ("maxR", maxR)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            # xNew = self.x - XCOM2
            # yNew = self.y - XCOM2
            # zNew = self.z - XCOM2
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units usint astropy and round all values
        COMP = np.around(COMP,2)
        COMP *= u.kpc

        # and then return the COM positon vector
        return COMP
        
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x - COMX.value
        yV = self.y - COMY.value
        zV = self.z - COMZ.value
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV <= RVMAX.value)

        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew  = self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew,mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # write your own code below
        COMV = [VXCOM, VYCOM, VZCOM]
        COMV = np.around(COMV,2)
        COMV = COMV * u.km / u.s

        # return the COM vector                                                                                        
        return COMV
    


# # ANSWERING QUESTIONS
# #######################


# # Create  a Center of mass object for the MW, M31 and M33
# # below is an example of using the class for MW
# MWCOM = CenterOfMass("MW_000.txt", 2)

# # below gives you an example of calling the class's functions
# # MW:   store the position and velocity COM 
# MW_COMP = MWCOM.COM_P(0.1)
# MW_COMV = MWCOM.COM_V(MW_COMP[0],MW_COMP[1],MW_COMP[2])

# # now write your own code to answer questions

# # Create objects for the other galaxies, extract their COM data
# M31COM = CenterOfMass("M31_000.txt", 2)
# M31_COMP = M31COM.COM_P(0.1)
# M31_COMV = M31COM.COM_V(M31_COMP[0],M31_COMP[1],M31_COMP[2])

# M33COM = CenterOfMass("M33_000.txt", 2)
# M33_COMP = M33COM.COM_P(0.1)
# M33_COMV = M33COM.COM_V(M33_COMP[0],M33_COMP[1],M33_COMP[2])


# # print COM data
# print("MW  COMP = {}".format(MW_COMP))
# print("MW  COMV = {}".format(MW_COMV))

# print("M31 COMP = {}".format(M31_COMP))
# print("M31 COMV = {}".format(M31_COMV))

# print("M33 COMP = {}".format(M33_COMP))
# print("M33 COMV = {}".format(M33_COMV))

# # calulate separation between galaxies
# sep1 = np.absolute(MW_COMP - M31_COMP)
# sep2 = np.absolute(M33_COMP - M31_COMP)

# # calulate separation between galaxies
# vel1 = np.absolute(MW_COMV - M31_COMV)
# vel2 = np.absolute(M33_COMV - M31_COMV)

# # print separations
# print("Magnitude of Separation between MW and M31:  {}".format(np.around(np.linalg.norm(sep1),2)))
# print("Magnitude of Separation between M33 and M31: {}".format(np.around(np.linalg.norm(sep2),2)))

# print("Magnitude of velocity between MW and M31:  {}".format(np.around(np.linalg.norm(vel1),2)))
# print("Magnitude of velocity between M33 and M31: {}".format(np.around(np.linalg.norm(vel2),2)))

# # print answer to Question 4
# print("\nQuestion 4: The iterative process used to determine the COM is important because as the galaxies merge, they begin to significantly affect each other and transfer material (this is the definition of merging). An iterative process that we have written continually lets us calculate the center of mass of all the particles within a radius, so we can know the COM infomation even though the galaxies may have merged.")














