# test to see if np.save() worked

# try and load a data from a .npy file
# load module
import numpy as np

# import data
data = np.load("LowResNPY/M31_004.npy")

# print to see if it worked
print(data)