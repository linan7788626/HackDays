#/usr/bin/env python
#from pylab import *

import numpy as np
import sys
#---------------------------------------------------------------------------------
import ctypes as ct
sps = ct.CDLL("./lib/libsphsdens.so")
sps.cal_sph_sdens.argtypes =[ct.c_char_p,ct.c_float,ct.c_long, \
                            ct.c_float,ct.c_long,ct.c_long, \
                            ct.c_float,ct.c_float,ct.c_float,ct.c_float, \
                            np.ctypeslib.ndpointer(dtype = ct.c_float), \
                            np.ctypeslib.ndpointer(dtype = ct.c_float), \
                            np.ctypeslib.ndpointer(dtype = ct.c_float)]

sps.cal_sph_sdens.restype  = ct.c_int

def call_sph_sdens(file_particles,Bsz,Nc,Np):

    dsx = ct.c_float(Bsz/Nc)
    Ngb = ct.c_long(32)
    xc1 = ct.c_float(0.0)
    xc2 = ct.c_float(0.0)
    xc3 = ct.c_float(0.0)
    mass_particle = ct.c_float(1.0)
    posx1 = np.zeros((Nc,Nc),dtype=ct.c_float)
    posx2 = np.zeros((Nc,Nc),dtype=ct.c_float)
    sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    sps.cal_sph_sdens(file_particles,ct.c_float(Bsz),ct.c_long(Nc),dsx,Ngb,ct.c_long(Np),xc1,xc2,xc3,mass_particle,posx1,posx2,sdens)
    return posx1,posx2,sdens
#---------------------------------------------------------------------------------
def make_r_coor(nc,dsx):

    bsz = nc*dsx
    x1 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0

    x2,x1 = np.meshgrid(x1,x2)
    return x1,x2
def make_c_coor(nc,dsx):

    bsz = nc*dsx
    x1,x2 = np.mgrid[0:(bsz-dsx):nc*1j,0:(bsz-dsx):nc*1j]-bsz/2.0+dsx/2.0
    return x1,x2

def ksz_based_on_sdens(zbar,vbar,sdens):
    f_b = 0.046/0.31
    sigma_T = 1.0
    vc = 3e5 #km/s
    a_scale = 1.0/(1.0+zbar)
    x_H = 1.0-y_He
    m_H = 1.0
    y_He = 1.0
    m_He = 1.0
    mu_xy = x_H/m_H+y_He/m_He
    res = -f_b*sigma_T*mu_xy/vc/a_scale*vbar*sdens
    return res

def tsz_based_on_sdens(zbar,vbar,sdens):
    f_b = 0.046/0.31
    sigma_T = 1.0
    vc = 3e5 #km/s
    a_scale = 1.0/(1.0+zbar)
    x_H = 1.0
    m_H = 1.0
    y_He = 1.0
    m_He = 1.0
    mu_xy = x_H/m_H+y_He/m_He
    res = -f_b*sigma_T*mu_xy/vc/a_scale*vbar*sdens
    return res
#---------------------------------------------------------------------------------
#from mpi4py import MPI
#@profile
def main():
    zbar = 0.7
    vbar = 500 #km/s
    ncc=1024
    npp=200000
    bsz=3.0

    file_particles = sys.argv[1]

    posx1,posx2,sdens = call_sph_sdens(file_particles,bsz,ncc,npp)
    #ksz_map = ksz_based_on_sdens(zbar,vbar,sdens)
    #tsz_map = tsz_based_on_sdens(zbar,vbar,sdens)

    import pylab as pl

    a0 = np.fromfile("../output/cnfw_ana_dense.bin",dtype=np.float32)
    a0 = a0.reshape((ncc,ncc))

    pl.figure()
    levels = [3.,4.,5.,6.,7.]

    pl.contour(posx1,posx2,np.log10(sdens),levels,colors=('r',))
    pl.contour(posx2,posx1,np.log10(a0),levels,colors=('k',))

    pl.colorbar()
    pl.show()

    return 0
main()
