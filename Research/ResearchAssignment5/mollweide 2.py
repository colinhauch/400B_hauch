import numpy as np
import matplotlib.pyplot as plt
import ephem # to make coordinate systems conversions
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
    ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
    ax.scatter(np.radians(x),np.radians(Dec))  # convert degrees to radians
    ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("RA")
    ax.xaxis.label.set_fontsize(12)
    ax.set_ylabel("Dec")
    ax.yaxis.label.set_fontsize(12)
    ax.grid(True)
np.random.seed(13)    # To make reproducible the plot
nrand = np.random.rand(100,2)    # Array (100,2) of random values in [0,1)
nrand *= np.array([360.,180.]) # array in [0,360) x[0,180) range
nrand -= np.array([0,90.]) # array in [0,360)x[-90,90) range
nrand = nrand[(-86 < nrand[:,1]) & (nrand[:,1] < 86)] # To avoid Matplotlib Runtime Warning,
RA = nrand[:,0]
Dec = nrand[:,1]
plot_mwd(RA,Dec)
plt.show()