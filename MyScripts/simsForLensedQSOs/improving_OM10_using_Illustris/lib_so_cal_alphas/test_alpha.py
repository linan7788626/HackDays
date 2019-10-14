#/usr/bin/env python
#from pylab import *

import numpy as np
import pyfits
import sys
#---------------------------------------------------------------------------------
import ctypes as ct
gls = ct.CDLL("./libglsg.so")
gls.sdens_to_alphas.argtypes = [ct.c_double,\
								np.ctypeslib.ndpointer(dtype = ct.c_double), \
								ct.c_int,ct.c_double,ct.c_double,ct.c_double,\
								np.ctypeslib.ndpointer(dtype = ct.c_double), \
								np.ctypeslib.ndpointer(dtype = ct.c_double)]
gls.sdens_to_alphas.restype  = ct.c_void_p

def call_cal_alphas(pmass,sdens,Ncc,boxsize,zl,zs):

	alpha1 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)
	alpha2 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)

	gls.sdens_to_alphas(pmass,sdens,Ncc,boxsize,zl,zs,alpha1,alpha2)

	return alpha1,alpha2

#---------------------------------------------------------------------------------
bls = ct.CDLL("./libbuildlos.so")
bls.make_los_images.argtypes = [ct.c_char_p,ct.c_int,\
								np.ctypeslib.ndpointer(dtype = ct.c_double), \
								np.ctypeslib.ndpointer(dtype = ct.c_double), \
								np.ctypeslib.ndpointer(dtype = ct.c_double)]
bls.make_los_images.restype  = ct.c_void_p

def call_build_los(filename,Ncc):

	ti = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)
	tv = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)
	tb = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)

	bls.make_los_images(filename,Ncc,ti,tv,tb)
	return ti,tv,tb

#---------------------------------------------------------------------------------
from mpi4py import MPI
@profile
def main():

	#nfiner = 4.0
	NCC=2048
	#NCC=4096
	DSL=0.09 #arcsec
	BSZ=DSL*NCC
	ZZL=0.55
	ZS0=10.0

	num_halos = int(sys.argv[1])
	sdens = np.fromfile(str(sys.argv[2]),dtype=np.float64)
	sdens = sdens.reshape((NCC,NCC))


	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()

	for i in xrange(rank,num_halos,size):
		call_cal_alphas(1.0,sdens,NCC,BSZ,ZZL,ZS0)

		#ti,tv,tb = call_build_los("./gals_in_HUDF/fits.list",NCC)

		#file_out = "los_images/bdded_"+str(i)+"_I.fits"
		#pyfits.writeto(file_out,ti,clobber=True)
		#file_out = "los_images/bdded_"+str(i)+"_V.fits"
		#pyfits.writeto(file_out,tv,clobber=True)
		#file_out = "los_images/bdded_"+str(i)+"_B.fits"
		#pyfits.writeto(file_out,tb,clobber=True)

	return 0
main()
