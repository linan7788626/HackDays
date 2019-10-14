#!/usr/bin/env python
import numpy as np
import ctypes as ct
#-------------------------------------------------------------
llt = ct.CDLL("./libcl_test.so")
#llt.main_test.argtypes = []
#llt.main_test.restype  = ct.c_int

llt.opencl_lq_py.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_float), \
                             np.ctypeslib.ndpointer(dtype =  ct.c_float), \
                             np.ctypeslib.ndpointer(dtype =  ct.c_float), \
                             ct.c_int,\
                             np.ctypeslib.ndpointer(dtype = ct.c_float), \
                             np.ctypeslib.ndpointer(dtype = ct.c_float), \
                             np.ctypeslib.ndpointer(dtype = ct.c_float), \
                             np.ctypeslib.ndpointer(dtype = ct.c_float)]
llt.opencl_lq_py.restype  = ct.c_int

#def call_llt_main():
#    err = ct.c_int(0)
#    err = llt.main_test()
#    return err

def call_opencl_lq(posx1,posx2,lpar):
    posx1 = np.array(posx1,dtype=ct.c_float)
    posx2 = np.array(posx2,dtype=ct.c_float)

    lpar = np.array(lpar,dtype=ct.c_float)
    count = len(posx1.flatten())

    alpha1 = np.zeros((count),dtype=ct.c_float)
    alpha2 = np.zeros((count),dtype=ct.c_float)
    alpha1_c = np.zeros((count),dtype=ct.c_float)
    alpha2_c = np.zeros((count),dtype=ct.c_float)

    llt.opencl_lq_py(posx1,posx2,lpar,count,alpha1,alpha2,alpha1_c,alpha2_c)
    return alpha1,alpha2,alpha1_c,alpha2_c

def main():
    xlc0 = 0.0
    ylc0 = 0.0
    ql0 = 0.7
    rc0 = 0.1
    re0 = 1.0
    phi0 = 0.0
    lpar = np.array([ylc0,xlc0,ql0,rc0,re0,phi0])

    count = 1024*1024*32
    xi1 = np.random.random(count)
    xi2 = np.random.random(count)
    alpha1,alpha2,alpha1_c,alpha2_c = call_opencl_lq(xi1,xi2,lpar)
    print alpha1,alpha2,alpha1_c,alpha2_c
    return 0

if __name__ == '__main__':
    main()
