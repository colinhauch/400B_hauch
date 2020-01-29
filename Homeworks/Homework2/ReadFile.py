# Colin Hauch
# 400B, with Dr. Gurtina Besla
# Homework2


# Location:
# /Users/chauch/Documents/College/Undergraduate/Semester8/400B/400B_hauch/Homeworks/Homework2

# import relavent modules
import numpy as np
import astropy.units as u

# Style choice: all variables are camelcased
# global variables are capitalized
# local variables are lowercased

# Define file name
# this is a hardcoded file name -- I'm not a fan
# could make the file name an input in the command to run this script
Filename = "MW_000.txt"

# define the function Read that takes the name of the file as an input
def Read(filename):
	# open the file
	file = open(filename, 'r')

	# read the first two lines of the file
	# "file.readline()" reads sequential lines,
	# Therefore, (line1 != line2)
	line1 = file.readline()
	line2 = file.readline()

	# split lines by white space
	label1, value1 = line1.split()
	label2, value2 = line2.split()

	# create values with astropy units
	time = float(value1)*u.Myr
	total = float(value2) # unit? unitless quantity?

	# # Test print statements for development
	# print("time: " + str(time))
	# print("total particles: " + str(total))

	# close the file: best practices
	file.close()

	# store the remainder of the file in this particular format
	# 
	# info from homework:
	# parameters: “dtype=None” means line is split using white spaces
	# “skip header=3” skipping the first 3 lines
	# the flag “names=True” creates arrays to store the data with the right labels 
	# the labels are : ”type, m, x, y,z, vx, vy, vz”
	data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)

	# return values and data for future use
	return time, total, data

	# label, value = line1.split()
	# time = float(value)*10.0*u.Myr

# collect united values from header of the file
Time, Total, Data = Read(Filename)


# # test print of certain data values:
# # For example, if you wanted to know the particle type of the 2nd particle in the file you would write:
# #
# # print(Data[’type’][1])
# #
# print(Data['type'][1])









