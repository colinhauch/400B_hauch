

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
#import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass import CenterOfMass




def OrbitCOM(galaxy, start, end, n):
    # This function will compute the time and COM position and velocity 
    # vectors of a given galaxy in each snapshot and save that output into a file
    #
    # INPUTS:
    #   galaxy - the name of the galaxy as a string, e.g. “MW”
    #   start - The number of the first snapshot to be read in   
    #   end - the number of the last snapshot to be read in.
    #   n - n integer indicating the intervals over which COM will be returned.
    #
    # RETURNS:
    #   No returns?
    
    # compose the filename for output
    fileout = "Orbit_"+galaxy+".txt"
    
    # set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    if galaxy == "M33":
        delta = 0.1
        VolDec = 4.0
    else:
        delta = 0.1
        VolDec = 2.0

    
    # generate the snapshot id sequence 
    snap_ids = np.arange(start,end,n)

    # it is always a good idea to also check if the input is eligible (not required) 
    if len(snap_ids) == 0:
        print("snap_ids array initialization failed. Check inputs")
        return
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    # NOTE: maybe could replace "len(snap_ids)" with "n", but not certain
    orbit = np.zeros((len(snap_ids),7))
    #orbit = np.zeros((n,7))
    
    # print for awareness of what is going on
    print("\n\ncalculating data for "+galaxy+"\n")

    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        loc = "/Users/colinhauch/Documents/College/Undergraduate/Semester8/400B/400B_hauch/Homeworks/Homework6"
        galaxyFolder = "/"+galaxy+"_VLowRes/"

        # create snap number string
        ilbl = '000' + str(snap_id)

        #remove all but last 3 digits to finalize snap number string
        ilbl = ilbl[-3:]

        # finalize the file name
        filename = loc+galaxyFolder+galaxy+"_"+ilbl+".txt"


        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename,2)

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        GalCOMP = COM.COM_P(delta,VolDec)
        #GalCOMV = COM.COM_V(GalCOMP[0].value, GalCOMP[1].value, GalCOMP[2].value)
        GalCOMV = COM.COM_V(GalCOMP[0], GalCOMP[1], GalCOMP[2])
    
    
        # store the time, pos, vel in ith element of the orbit array, without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i] = COM.time.value/1000, GalCOMP[0].value, GalCOMP[1].value, GalCOMP[2].value, GalCOMV[0].value, GalCOMV[1].value, GalCOMV[2].value

        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 


# OrbitCOM("MW",0,800,5)
# OrbitCOM("M31",0,800,5)
# OrbitCOM("M33",0,800,5)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt

MW_data = np.genfromtxt("Orbit_MW.txt",dtype=None,names=True,skip_header=0)
M31_data = np.genfromtxt("Orbit_M31.txt",dtype=None,names=True,skip_header=0)
M33_data = np.genfromtxt("Orbit_M33.txt",dtype=None,names=True,skip_header=0)



# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def vectorDiff(v1,v2):
    # Calculates v1 - v2 
    # INPUTS:
    #   v1 - first vector
    #   v2 - second vector
    #
    # RETURNS:
    #   vector difference between input vectors
    return np.sqrt((v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2)


# Determine the magnitude of the relative position and velocities 

# of MW and M31
MW_M31_POS = vectorDiff((MW_data['x'],MW_data['y'],MW_data['z']),(M31_data['x']-MW_data['x'],M31_data['y']-MW_data['y'],M31_data['z']-MW_data['z']))
MW_M31_VEL = vectorDiff((MW_data['vx'],MW_data['vy'],MW_data['vz']),(M31_data['vx']-MW_data['vx'],M31_data['vy']-MW_data['vy'],M31_data['vz']-MW_data['vz']))
# of M33 and M31
M33_M31_POS = vectorDiff((M33_data['x'],M33_data['y'],M33_data['z']),(M31_data['x']-M33_data['x'],M31_data['y']-M33_data['y'],M31_data['z']-M33_data['z']))
M33_M31_VEL = vectorDiff((M33_data['vx'],M33_data['vy'],M33_data['vz']),(M31_data['vx']-M33_data['vx'],M31_data['vy']-M33_data['vy'],M31_data['vz']-M33_data['vz']))




# Plot the Orbit of the galaxies 
#################################

fig = plt.figure(figsize = (10,10))
ax = plt.subplot(111)

#Plot mag diff as a function of time (column 0)
plt.plot(MW_data['t'], MW_M31_POS,color='blue',linewidth=5,label='MW_M31')

#axis labels
plt.xlabel('Time (Gyr)', fontsize = 22 )
plt.ylabel('Separation (kpc)', fontsize = 22)

#add legend
#legend = ax.legend(loc = 'lower right', fontsize = 'x-large')
plt.legend()
plt.show()

fig = plt.figure(figsize = (10,10))
ax = plt.subplot(111)

#Plot mag diff as a function of time (column 0)
plt.plot(M33_data['t'], M33_M31_POS,color='blue',linewidth=5,label='M33_M31')

#axis labels
plt.xlabel('Time (Gyr)', fontsize = 22 )
plt.ylabel('Separation (kpc)', fontsize = 22)

#add legend
#legend = ax.legend(loc = 'lower right', fontsize = 'x-large')
plt.legend()
plt.show()


# Plot the orbital velocities of the galaxies 
#################################
fig = plt.figure(figsize = (10,10))
ax = plt.subplot(111)

#Plot mag diff as a function of time (column 0)
plt.plot(MW_data['t'], MW_M31_VEL,color='red',linewidth=5,label='MW_M31')

#axis labels
plt.xlabel('Time (Gyr)', fontsize = 22 )
plt.ylabel('relative velocity (km/s)', fontsize = 22)

#add legend
#legend = ax.legend(loc = 'lower right', fontsize = 'x-large')
plt.legend()
plt.show()


fig = plt.figure(figsize = (10,10))
ax = plt.subplot(111)

#Plot mag diff as a function of time (column 0)
plt.plot(MW_data['t'], M33_M31_VEL,color='red',linewidth=5,label='M33_M31')

#axis labels
plt.xlabel('Time (Gyr)', fontsize = 22 )
plt.ylabel('relative velocity (km/s)', fontsize = 22)

#add legend
#legend = ax.legend(loc = 'lower right', fontsize = 'x-large')
plt.legend()
plt.show()


"""
QUESTIONS:

1. MW and M31 will have 2 close encounters before merging.

2. When the galaxies of interest are further apart, their relative velocities
    are smaller. When they approach each other, their relative velocities
    increase -- this effect causes them to shoot past / through one another.

3. M31 and MW merge at t=~6 Gyr. M33's orbit decays when they merge.
"""


