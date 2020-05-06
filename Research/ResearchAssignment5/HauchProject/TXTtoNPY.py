# Colin Hauch
# This script creates ".npy" files for each text datafile
# so that they don't have to be read every time the a 
# subsequent file is run

# import relevant modules
from ReadFile import Read
import numpy as np

# initialize source file parameters
sourceFolder = "LowRes/"
sourceSuffix = ".txt"

productFolder = "LowResNPY/"
productSuffix = ".npy"

# snap number from 000 to 801

# loop through each filename
for i in range(802):

	# gen snap number from 000 to 801
	if len(str(i)) == 3:
		snap = str(i)
	elif len(str(i)) == 2:
		snap = "0" + str(i)
	else:
		snap = "00" + str(i)

	# construct source filenames
	sMWfile  = sourceFolder+"MW_"+snap+sourceSuffix
	sM31file = sourceFolder+"M31_"+snap+sourceSuffix
	sM33file = sourceFolder+"M33_"+snap+sourceSuffix

	pMWfile   = productFolder+"MW_"+snap+productSuffix
	pM31file  = productFolder+"M31_"+snap+productSuffix
	pM33file  = productFolder+"M33_"+snap+productSuffix


	# read in text files
	MWtime, MWtotal, MWdata 	= Read(sMWfile)
	M31time, M31total, M31data 	= Read(sM31file)
	M33time, M33total, M33data 	= Read(sM33file)


	# save npy files
	np.save(pMWfile,MWdata)
	np.save(pM31file,M31data)
	np.save(pM33file,M33data)

	# print snap for progress updates
	print(snap)



