#!/usr/bin/env python

import numpy as np
import ctypes as ct

Dtype=ct.c_double
#Dtype=np.float32_t

#-------------------------------------------------------------

cv_test = ct.CDLL("./liball_cv_test.so")
cv_test.lens_images.argtypes =[np.ctypeslib.ndpointer(dtype =  Dtype),\
                               np.ctypeslib.ndpointer(dtype =  Dtype), \
                               ct.c_int,ct.c_int,\
                               np.ctypeslib.ndpointer(dtype = Dtype),ct.c_int,\
                               np.ctypeslib.ndpointer(dtype = Dtype),ct.c_int,\
                               np.ctypeslib.ndpointer(dtype = Dtype)]
cv_test.lens_images.restype  = ct.c_void_p

def call_lens_images(xi1,xi2,gpar,gpars):
    nx1,nx2 = np.shape(xi1)
    npars = len(gpar)
    nsubs = len(gpars)

    xi1 = np.array(xi1,dtype=Dtype)
    xi2 = np.array(xi2,dtype=Dtype)
    gpar_array = np.array(gpar,dtype=Dtype)
    gpars_array = np.array(gpars,dtype=Dtype)
    g_lens = np.zeros((nx1,nx2),dtype=Dtype)

    cv_test.lens_images(xi1,xi2,ct.c_int(nx1),ct.c_int(nx2),gpar_array,ct.c_int(npars),gpars_array,ct.c_int(nsubs),g_lens)
    return g_lens

#-------------------------------------------------------------
cv_test.mmbr_images.argtypes =[np.ctypeslib.ndpointer(dtype =  Dtype),\
                               np.ctypeslib.ndpointer(dtype =  Dtype), \
                               ct.c_int,ct.c_int,\
                               np.ctypeslib.ndpointer(dtype = Dtype),ct.c_int,\
                               np.ctypeslib.ndpointer(dtype = Dtype),ct.c_int,\
                               np.ctypeslib.ndpointer(dtype = Dtype)]
cv_test.mmbr_images.restype  = ct.c_void_p

def call_mmbr_images(xi1,xi2,gpar,gpars):
    nx1,nx2 = np.shape(xi1)
    npars = len(gpar)
    nsubs = len(gpars)

    xi1 = np.array(xi1,dtype=Dtype)
    xi2 = np.array(xi2,dtype=Dtype)
    gpar_array = np.array(gpar,dtype=Dtype)
    gpars_array = np.array(gpars,dtype=Dtype)
    g_edge = np.zeros((nx1,nx2),dtype=Dtype)

    cv_test.mmbr_images(xi1,xi2,ct.c_int(nx1),ct.c_int(nx2),gpar_array,ct.c_int(npars),gpars_array,ct.c_int(nsubs),g_edge)

    return g_edge
#-------------------------------------------------------------
cv_test.all_about_lensing.argtypes =[np.ctypeslib.ndpointer(dtype =  Dtype),\
                                     np.ctypeslib.ndpointer(dtype =  Dtype), \
                                     ct.c_int,ct.c_int,\
                                     np.ctypeslib.ndpointer(dtype = Dtype),\
                                     np.ctypeslib.ndpointer(dtype = Dtype),ct.c_int,\
                                     np.ctypeslib.ndpointer(dtype = Dtype),ct.c_int,\
                                     np.ctypeslib.ndpointer(dtype = Dtype),\
                                     np.ctypeslib.ndpointer(dtype = Dtype),\
                                     np.ctypeslib.ndpointer(dtype = Dtype),\
                                     np.ctypeslib.ndpointer(dtype = Dtype)]
cv_test.all_about_lensing.restype  = ct.c_void_p


def call_all_about_lensing(xi1,xi2,spar,lpar,lpars):

    nx1,nx2 = np.shape(xi1)
    npars = len(lpar)
    nsubs = len(lpars)

    xi1 = np.array(xi1,dtype=Dtype)
    xi2 = np.array(xi2,dtype=Dtype)
    spar = np.array(spar,dtype=Dtype)
    lpar = np.array(lpar,dtype=Dtype)
    lpars = np.array(lpars,dtype=Dtype)

    s_image = np.zeros((nx1,nx2),dtype=Dtype)
    g_lensimage = np.zeros((nx1,nx2),dtype=Dtype)
    critical = np.zeros((nx1,nx2),dtype=Dtype)
    caustic = np.zeros((nx1,nx2),dtype=Dtype)

    cv_test.all_about_lensing(xi1,xi2,ct.c_int(nx1),ct.c_int(nx2),spar,lpar,ct.c_int(npars),lpars,ct.c_int(nsubs),s_image,g_lensimage,critical,caustic)

    return s_image,g_lensimage,critical,caustic.T
