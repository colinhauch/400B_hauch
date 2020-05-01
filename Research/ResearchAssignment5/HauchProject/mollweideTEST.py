import numpy as np
import matplotlib.pyplot as plt
import ephem # to make coordinate systems conversions

# load list of INDEXES of particles that could be used as VPP
cannidateParticles = np.load("cannidates.npy")

# pick one of the indexes
# I picked the 30th particle, as a proof of concept
vppIndex = cannidateParticles[30]

#print(vppIndex)

# load snap data
snap000Data = np.load("LowResNPY/M31_000.npy")

# filter for disk particles (type 2)
index = np.where(snap000Data['type'] == 2)

# fetch the line of data corresponding to VPP
vpp = snap000Data[vppIndex]

# print data to make sure I'm not crazy and it works
print(vpp)

# fetch all data
# should possibly remove VPP particle...
# but that would change the index of all
# the particles after and numpy seems to 
# handle the divByZero error within the
# trig functions
xdata = snap000Data['x'][index] - vpp[2]
ydata = snap000Data['y'][index] - vpp[3]
zdata = snap000Data['z'][index] - vpp[4]

# latitude = np.arctan(zdata/np.sqrt(xdata**2 + ydata**2))
# #longitude = np.arctan(ydata/xdata)

# why convert lat,long to RA,Dec when I can calculate them directly?
RA = np.arctan2(ydata,xdata) * 180 / np.pi
DEC = np.arcsin(zdata/np.sqrt(xdata**2 + ydata**2 + zdata**2)) * 180 / np.pi

# test print
print(RA)
print(DEC)


# mollweid plotting method from 
# http://balbuceosastropy.blogspot.com/2013/09/the-mollweide-projection.html
def plot_mwd(RA,Dec,org=0,title='Mollweide projection', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    x = np.remainder(RA+360-org,360) # shift RA values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+org,360)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111, projection=projection)
    #ax.scatter(np.radians(x),np.radians(Dec))  # convert degrees to radians
    ax.hist2d(np.radians(x),np.radians(Dec))
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
# np.random.seed(13)    # To make reproducible the plot
# nrand = np.random.rand(100,2)    # Array (100,2) of random values in [0,1)
# nrand *= np.array([360.,180.]) # array in [0,360) x[0,180) range
# nrand -= np.array([0,90.]) # array in [0,360)x[-90,90) range
# nrand = nrand[(-86 < nrand[:,1]) & (nrand[:,1] < 86)] # To avoid Matplotlib Runtime Warning,
# RA = nrand[:,0]
# Dec = nrand[:,1]

# plot data
plot_mwd(RA,DEC)

# show plot!
plt.show()