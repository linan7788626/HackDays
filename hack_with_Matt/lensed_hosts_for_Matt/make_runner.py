import numpy as np
import sys

the_directory="/Users/benedictine/Documents/lsst/instance_files/"


b=np.loadtxt(the_directory+"instance.list", dtype="str")

name = np.loadtxt(
        the_directory+"instance.list", dtype="str", comments='#', 
        converters=None, skiprows=0, usecols=(0),
        unpack=True, ndmin=0)

phile = open(the_directory+"runner.sh","w")

for p in range (0, len(b)):
	print >> phile, "./phosim", "instance_files/"+name[p], "-c trim.txt -w workdir/work"+str(p), "-o outdir/out"+str(p) 

phile.close()