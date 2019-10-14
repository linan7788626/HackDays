#/usr/bin/env python
import sys
del sys.path[2]
import numpy as np
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
    mass_particle = ct.c_float(1e9)
    posx1 = np.zeros((Nc,Nc),dtype=ct.c_float)
    posx2 = np.zeros((Nc,Nc),dtype=ct.c_float)
    sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    sps.cal_sph_sdens(file_particles,ct.c_float(Bsz),ct.c_long(Nc),dsx,Ngb,ct.c_long(Np),xc1,xc2,xc3,mass_particle,posx1,posx2,sdens)
    return sdens

def ksz_based_on_sdens(zbar,vbar,sdens):
    f_b = 0.0486/0.3089
    sigma_T = 0.665246e-24 # cm^2
    vc = 3e5 #km/s
    a_scale = 1.0/(1.0+zbar)
    y_He = 0.2477
    x_H = 1.0-y_He
    m_H = 1.67372e-27 #kg
    m_He = 6.64647e-27 #k g
    #T_cmb = 2.7255 # k
    mu_xy = x_H/m_H+y_He/m_He
    norm_factor = 1.5e-19 # Msun/h/kg * (cm/Mpc/h)^2
    res = -f_b*sigma_T*mu_xy/vc/a_scale*vbar*sdens*norm_factor
    return res #Delta_T/T_cmb

def tsz_based_on_sdens(zbar,vbar,sdens):
    '''
    Need to ask Lindsey.
    '''
    return
#---------------------------------------------------------------------------------
#from mpi4py import MPI
#@profile
def main():
    zbar = 0.7
    vbar = 500 #km/s comoving distance
    ncc=1024
    npp=200000 # total mass = 1e9*2e5 M_sun/h
    bsz=3.0 #Mpc/h comoving distance

    file_particles = sys.argv[1]

    sdens = call_sph_sdens(file_particles,bsz,ncc,npp)
    sdens.astype(np.float32).tofile("./output_files/cnfw_sdens.bin")

    ksz_map = ksz_based_on_sdens(zbar,vbar,sdens)
    #---------------------------
    # Save 2d array to a binary file.
    ksz_map.astype(np.float32).tofile("./output_files/cnfw_ksz.bin")
    #---------------------------


    #---------------------------
    # Read the output binary files.
    a0 = np.fromfile("./output_files/cnfw_ksz.bin",dtype=np.float32)
    a0 = a0.reshape((ncc,ncc))
    #---------------------------


    #---------------------------
    # Plot the contours
    import pylab as pl
    pl.contourf(a0)
    pl.colorbar()
    pl.show()


    return 0
main()
