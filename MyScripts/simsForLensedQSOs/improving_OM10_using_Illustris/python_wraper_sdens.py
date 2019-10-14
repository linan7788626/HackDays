#/usr/bin/env python
import sys
del sys.path[2]
import numpy as np
import pylab as pl
import illustris_python as il
basePath = '../Illustris-3'
#---------------------------------------------------------------------------------
import ctypes as ct
sps = ct.CDLL("../calsdens_sph_libs/lib/libsphsdens.so")
#sps.cal_sph_sdens.argtypes =[ct.c_char_p,ct.c_float,ct.c_long, \
                            #ct.c_float,ct.c_long,ct.c_long, \
                            #ct.c_float,ct.c_float,ct.c_float,ct.c_float, \
                            #np.ctypeslib.ndpointer(dtype = ct.c_float), \
                            #np.ctypeslib.ndpointer(dtype = ct.c_float), \
                            #np.ctypeslib.ndpointer(dtype = ct.c_float)]

#sps.cal_sph_sdens.restype  = ct.c_int

#def call_sph_sdens(file_particles,Bsz,Nc,Np,pmass):

    #dsx = ct.c_float(Bsz/Nc)
    #Ngb = ct.c_long(32)
    #xc1 = ct.c_float(0.0)
    #xc2 = ct.c_float(0.0)
    #xc3 = ct.c_float(0.0)
    #mass_particle = ct.c_float(pmass)
    #posx1 = np.zeros((Nc,Nc),dtype=ct.c_float)
    #posx2 = np.zeros((Nc,Nc),dtype=ct.c_float)
    #sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    #sps.cal_sph_sdens(file_particles,ct.c_float(Bsz),ct.c_long(Nc),dsx,Ngb,ct.c_long(Np),xc1,xc2,xc3,mass_particle,posx1,posx2,sdens)
    #return sdens

sps.cal_sph_sdens_arrays.argtypes =[np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     ct.c_float,ct.c_long, \
                                     ct.c_float,ct.c_long,ct.c_long, \
                                     ct.c_float,ct.c_float,ct.c_float,ct.c_float, \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float)]

sps.cal_sph_sdens_arrays.restype  = ct.c_int

def call_sph_sdens_arrays(xi1,xi2,xi3,xic1,xic2,xic3,Bsz,Nc,Np,pmass):

    x1 = xi1.astype(ct.c_float)
    x2 = xi2.astype(ct.c_float)
    x3 = xi3.astype(ct.c_float)

    dsx = ct.c_float(Bsz/Nc)
    Ngb = ct.c_long(32)
    xc1 = ct.c_float(xic1)
    xc2 = ct.c_float(xic2)
    xc3 = ct.c_float(xic3)
    mass_particle = ct.c_float(pmass)
    posx1 = np.zeros((Nc,Nc),dtype=ct.c_float)
    posx2 = np.zeros((Nc,Nc),dtype=ct.c_float)
    sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    sps.cal_sph_sdens_arrays(x1,x2,x3,ct.c_float(Bsz),ct.c_long(Nc),dsx,Ngb,ct.c_long(Np),xc1,xc2,xc3,mass_particle,posx1,posx2,sdens)
    return sdens
#---------------------------------------------------------------------------------
if __name__ == '__main__':
    dm = il.snapshot.loadHalo(basePath,135,19,'dm')
    npp = dm['count']
    x1 = dm['Coordinates'][:,0]
    x2 = dm['Coordinates'][:,1]
    x3 = dm['Coordinates'][:,2]
    ncc=1024
    bsz=x1.max()-x1.min()
    print npp
    xc1 = np.mean(x1)
    xc2 = np.mean(x2)
    xc3 = np.mean(x3)
    pmass = 1e7
    sdens = call_sph_sdens_arrays(x1,x2,x3,xc1,xc2,xc3,bsz,ncc,npp,pmass)

    print np.shape(sdens)

