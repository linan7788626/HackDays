#!/usr/bin/env python
import numpy as np
import pylab as pl
import ctypes as ct

#--------------------------------------------------------------------
tri = ct.CDLL("../triangles_mapping/libtri.so")
tri.PIT.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
                    np.ctypeslib.ndpointer(dtype = ct.c_double), \
                    np.ctypeslib.ndpointer(dtype = ct.c_double), \
                    np.ctypeslib.ndpointer(dtype = ct.c_double)]
tri.PIT.restype  = ct.c_bool

def call_PIT(pt,v0,v1,v2):

    pt_in = np.array(pt,dtype=ct.c_double)
    v0_in = np.array(v0,dtype=ct.c_double)
    v1_in = np.array(v1,dtype=ct.c_double)
    v2_in = np.array(v2,dtype=ct.c_double)

    res = tri.PIT(pt_in,v0_in,v1_in,v2_in)
    return res

#--------------------------------------------------------------------
tri.Cart2Bary.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double)]
tri.Cart2Bary.restype  = ct.c_void_p

def call_cart_to_bary(pt,v0,v1,v2):

    pt_in = np.array(pt,dtype=ct.c_double)
    v0_in = np.array(v0,dtype=ct.c_double)
    v1_in = np.array(v1,dtype=ct.c_double)
    v2_in = np.array(v2,dtype=ct.c_double)

    bary_out = np.array([0,0,0],dtype=ct.c_double)
    tri.Cart2Bary(pt_in,v0_in,v1_in,v2_in,bary_out)
    return bary_out

#--------------------------------------------------------------------
tri.bary2cart.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double), \
                          np.ctypeslib.ndpointer(dtype = ct.c_double)]
tri.bary2cart.restype  = ct.c_void_p

def call_bary_to_cart(v0,v1,v2,bary):
    v0_in = np.array(v0,dtype=ct.c_double)
    v1_in = np.array(v1,dtype=ct.c_double)
    v2_in = np.array(v2,dtype=ct.c_double)
    bary_in = np.array(bary,dtype=ct.c_double)

    pt_out = np.array([0,0],dtype=ct.c_double)

    tri.bary2cart(v0_in,v1_in,v2_in,bary_in,pt_out)
    return pt_out

#--------------------------------------------------------------------
tri.mapping_triangles.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                  np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                  np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                  np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                  np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                  ct.c_int, \
                                  np.ctypeslib.ndpointer(dtype = ct.c_double)]
tri.mapping_triangles.restype  = ct.c_void_p

def call_mapping_triangles(pys,xi1,xi2,yi1,yi2):

    pys_in = np.array(pys,dtype=ct.c_double)
    xi1_in = np.array(xi1,dtype=ct.c_double)
    xi2_in = np.array(xi2,dtype=ct.c_double)
    yi1_in = np.array(yi1,dtype=ct.c_double)
    yi2_in = np.array(yi2,dtype=ct.c_double)
    nc_in = ct.c_int(np.shape(xi1)[0])

    xroots_out = np.zeros((10),dtype=ct.c_double)

    tri.mapping_triangles(pys_in,xi1_in,xi2_in,yi1_in,yi2_in,nc_in,xroots_out)

    return xroots_out
#--------------------------------------------------------------------
def tria_vectors(i,j,xgrids1,xgrids2):

    ip1 = i+1
    jp1 = j+1

    v1 = np.array([xgrids1[i,j],xgrids2[i,j]])
    v2 = np.array([xgrids1[ip1,j],xgrids2[ip1,j]])
    v3 = np.array([xgrids1[ip1,jp1],xgrids2[ip1,jp1]])
    return v1,v2,v3


def trib_vectors(i,j,xgrids1,xgrids2):

    ip1 = i+1
    jp1 = j+1

    v1 = np.array([xgrids1[i,j],xgrids2[i,j]])
    v2 = np.array([xgrids1[ip1,jp1],xgrids2[ip1,jp1]])
    v3 = np.array([xgrids1[i,jp1],xgrids2[i,jp1]])
    return v1,v2,v3


def mapping_triangles(ys1,ys2,xgrids1,xgrids2,lgrids1,lgrids2):

    ntris = np.shape(xgrids1)[0]-1
    xroots = []

    ncount = 0
    for i in xrange(ntris):
        for j in xrange(ntris):
            xv1,xv2,xv3 = tria_vectors(i,j,xgrids1,xgrids2)
            lv1,lv2,lv3 = tria_vectors(i,j,lgrids1,lgrids2)
            if call_PIT([ys1,ys2],lv1,lv2,lv3):
                baryc = call_cart_to_bary([ys1,ys2],lv1,lv2,lv3)
                xroots.append(call_bary_to_cart(xv1,xv2,xv3,baryc))
                ncount = ncount + 1

    for i in xrange(ntris):
        for j in xrange(ntris):
            xv1,xv2,xv3 = trib_vectors(i,j,xgrids1,xgrids2)
            lv1,lv2,lv3 = trib_vectors(i,j,lgrids1,lgrids2)
            if call_PIT([ys1,ys2],lv1,lv2,lv3):
                baryc = call_cart_to_bary([ys1,ys2],lv1,lv2,lv3)
                xroots.append(call_bary_to_cart(xv1,xv2,xv3,baryc))
                ncount = ncount + 1

    return np.array(xroots),ncount
#--------------------------------------------------------------------

def run_main():
    nnn = 512
    # bsz = 3600

    ys1 = 547.4
    ys2 = -1346.8

    xi1 = np.fromfile("./xi1.bin",dtype=np.float32).reshape((nnn,nnn))
    xi2 = np.fromfile("./xi2.bin",dtype=np.float32).reshape((nnn,nnn))
    yf1 = np.fromfile("./yf1.bin",dtype=np.float32).reshape((nnn,nnn))
    yf2 = np.fromfile("./yf2.bin",dtype=np.float32).reshape((nnn,nnn))

    # xroots = mapping_triangles(ys1,ys2,xi1,xi2,yf1,yf2)[0][0]
    # xroots = mapping_triangles(ys1,ys2,xi1,xi2,xi1,xi2)[0][0]
    # xroots = call_mapping_triangles([ys1,ys2],xi1,xi2,yf1,yf2)
    # xroots = call_mapping_triangles([ys1,ys2],xi1,xi2,xi1,xi2)

    print xroots

    # pl.figure(figsize=(10,10))
    # pl.plot(xroots[0],xroots[1],'go')
    # pl.plot(ys1,ys2,'ko')
    return 0

if __name__ == '__main__':
    run_main()
    pl.show()
