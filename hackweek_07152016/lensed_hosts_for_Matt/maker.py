import numpy as np
import sys

the_directory="/Users/benedictine/Documents/lsst/instance_files/"


b=np.loadtxt(the_directory+"instance.list", dtype="str")

name = np.loadtxt(
        the_directory+"instance.list", dtype="str", comments='#', 
        converters=None, skiprows=0, usecols=(0),
        unpack=True, ndmin=0)

phile1 = open(the_directory+"maker1.sh","w")
phile2 = open(the_directory+"maker2.sh","w")

for p in range (0, len(b)):
#	print >> phile, "./phosim", "instance_files/"+name[p], "-w work"+str(p), "-o out"+str(p) 
	print >> phile1, "mkdir out"+str(p)
	print >> phile2, "mkdir work"+str(p)

phile1.close()
phile2.close()