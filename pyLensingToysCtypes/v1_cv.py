#!/usr/bin/env python

import numpy as np
#import pylab as pl

import ctypes as ct
fcic = ct.CDLL("./libfcic.so")
fcic.forward_cic.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_double),\
                             np.ctypeslib.ndpointer(dtype =  ct.c_double), \
                             np.ctypeslib.ndpointer(dtype =  ct.c_double), \
                             ct.c_double,ct.c_double,ct.c_int,ct.c_int,ct.c_int,\
                             np.ctypeslib.ndpointer(dtype = ct.c_double)]
fcic.forward_cic.restype  = ct.c_void_p

def call_forward_cic(nx1,nx2,boxsize,yif1,yif2):
    img_in = np.array(np.ones(len(yif1)),dtype=ct.c_double)
    yif1 = np.array(yif1,dtype=ct.c_double)
    yif2 = np.array(yif2,dtype=ct.c_double)
    img_out = np.zeros((nx1,nx2))
    fcic.forward_cic(img_in,yif1,yif2,ct.c_double(boxsize),ct.c_double(boxsize),ct.c_int(nx1),ct.c_int(nx2),ct.c_int(len(yif1)),img_out)
    return img_out.T
#-------------------------------------------------------------

cv_test = ct.CDLL("./libcv_test.so")
cv_test.xy_rotate.argtypes =[np.ctypeslib.ndpointer(dtype =  ct.c_float),\
                             np.ctypeslib.ndpointer(dtype =  ct.c_float), \
                             ct.c_int,ct.c_int,ct.c_float,ct.c_float,ct.c_float,\
                             np.ctypeslib.ndpointer(dtype = ct.c_float),\
                             np.ctypeslib.ndpointer(dtype = ct.c_float)]
cv_test.xy_rotate.restype  = ct.c_void_p

def call_xy_rotate(x, y, xcen, ycen, phi):
    nx1,nx2 = np.shape(x)

    x = np.array(x,dtype=ct.c_float)
    y = np.array(y,dtype=ct.c_float)
    xnew = np.zeros((nx1,nx2),dtype=ct.c_float)
    ynew = np.zeros((nx1,nx2),dtype=ct.c_float)

    cv_test.xy_rotate(x,y,ct.c_int(nx1),ct.c_int(nx2),ct.c_float(xcen),ct.c_float(ycen),ct.c_float(phi),xnew,ynew)
    return (xnew,ynew)

#def call_xy_rotate(x, y, xcen, ycen, phi):
#
#    phirad = np.deg2rad(phi)
#    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
#    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
#    return (xnew,ynew)
##-------------------------------------------------------------
#cv_test.gauss_2d.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_float),\
#                             np.ctypeslib.ndpointer(dtype =  ct.c_float), \
#                             ct.c_int,ct.c_int,\
#                             np.ctypeslib.ndpointer(dtype = ct.c_float),\
#                             np.ctypeslib.ndpointer(dtype = ct.c_float)]
#cv_test.gauss_2d.restype  = ct.c_void_p
#
#def call_gauss_2d(x, y, par):
#    nx1,nx2 = np.shape(x)
#
#    x = np.array(x,dtype=ct.c_float)
#    y = np.array(y,dtype=ct.c_float)
#    par = np.array(par,dtype=ct.c_float)
#    res = np.zeros((nx1,nx2),dtype=ct.c_float)
#
#    cv_test.gauss_2d(x,y,ct.c_int(nx1),ct.c_int(nx2),par,res)
#
#    return res

def call_gauss_2d(x, y, par):
    #gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis])
    #gpars = np.asarray([aps,l_sigs,ylcs,xlcs,qls,phis])

    (xnew,ynew) = call_xy_rotate(x, y, par[0], par[1], par[5])
    r_ell_sq = ((xnew**2)*par[2] + (ynew**2)/par[2]) / np.abs(par[4])**2
    return par[3] * np.exp(-0.5*r_ell_sq)
##-------------------------------------------------------------
#cv_test.tophat_2d.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_float),\
#                             np.ctypeslib.ndpointer(dtype =  ct.c_float), \
#                             ct.c_int,ct.c_int,\
#                             np.ctypeslib.ndpointer(dtype = ct.c_float),\
#                             np.ctypeslib.ndpointer(dtype = ct.c_float)]
#cv_test.tophat_2d.restype  = ct.c_void_p
#
#def call_tophat_2d(x, y, par):
#
#    nx1,nx2 = np.shape(x)
#
#    x = np.array(x,dtype=ct.c_float)
#    y = np.array(y,dtype=ct.c_float)
#    par = np.array(par,dtype=ct.c_float)
#    res = np.zeros((nx1,nx2),dtype=ct.c_float)
#
#    cv_test.tophat_2d(x,y,ct.c_int(nx1),ct.c_int(nx2),par,res)
#    return res

def call_tophat_2d(x, y, par):
    (xnew,ynew) = call_xy_rotate(x, y, par[0], par[1], par[5])
    r_ell = np.sqrt(((xnew**2)*par[2] + (ynew**2)/par[2]))
    res = r_ell*0.0
    res[r_ell>=par[4]] = -1.0
    res[r_ell<par[4]] = 10000.0
    return res
#-------------------------------------------------------------
cv_test.lq_nie.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_float),\
                           np.ctypeslib.ndpointer(dtype =  ct.c_float), \
                           ct.c_int,ct.c_int,\
                           np.ctypeslib.ndpointer(dtype = ct.c_float),\
                           np.ctypeslib.ndpointer(dtype = ct.c_float),\
                           np.ctypeslib.ndpointer(dtype = ct.c_float)]
cv_test.lq_nie.restype  = ct.c_void_p

def call_lq_nie(x1,x2,lpar):
    nx1,nx2 = np.shape(x1)

    x1 = np.array(x1,dtype=ct.c_float)
    x2 = np.array(x2,dtype=ct.c_float)
    lpar = np.array(lpar,dtype=ct.c_float)
    res1 = np.zeros((nx1,nx2),dtype=ct.c_float)
    res2 = np.zeros((nx1,nx2),dtype=ct.c_float)

    cv_test.lq_nie(x1,x2,ct.c_int(nx1),ct.c_int(nx2),lpar,res1,res2)
    return res1,res2

#def call_lq_nie(x1,x2,lpar):
#
#    xc1 = lpar[0]
#    xc2 = lpar[1]
#    q   = lpar[2]
#    rc  = lpar[3]
#    re  = lpar[4]
#    pha = lpar[5]
#
#    phirad = np.deg2rad(pha)
#    cosa = np.cos(phirad)
#    sina = np.sin(phirad)
#
#    xt1 = (x1-xc1)*cosa+(x2-xc2)*sina
#    xt2 = (x2-xc2)*cosa-(x1-xc1)*sina
#
#    phi = np.sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc)
#    sq = np.sqrt(1.0-q*q)
#    pd1 = phi+rc/q
#    pd2 = phi+rc*q
#    fx1 = sq*xt1/pd1
#    fx2 = sq*xt2/pd2
#    qs = np.sqrt(q)
#
#    a1 = qs/sq*np.arctan(fx1)
#    a2 = qs/sq*np.arctanh(fx2)
#
#    res1 = (a1*cosa-a2*sina)*re
#    res2 = (a2*cosa+a1*sina)*re
#    return res1,res2
##--------------------------------------------------------------------
cv_test.find_critical_curve.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_float),\
                           ct.c_int,ct.c_int,\
                           np.ctypeslib.ndpointer(dtype = ct.c_float)]
cv_test.find_critical_curve.restype  = ct.c_void_p

def call_find_critical_curve(mu):
    nx1,nx2 = np.shape(mu)

    mu = np.array(mu,dtype=ct.c_float)
    res = np.zeros((nx1,nx2),dtype=ct.c_float)

    cv_test.find_critical_curve(mu,ct.c_int(nx1),ct.c_int(nx2),res)
    return res

#def call_find_critical_curve(mu):
#    #-------------------
#    # NEEDED
#    #-------------------
#    rows,cols = np.indices(np.shape(mu))
#    cdtn = np.sign(mu)*(np.sign(mu[rows-1,cols])+np.sign(mu[rows,cols-1])+np.sign(mu[(rows+1)%len(rows),cols])+np.sign(mu[rows,(cols+1)%len(cols)]))
#
#    res = mu*0.0
#    res[cdtn<4] = 1
#    res[cdtn>=4] = 0
#
#    return res
#--------------------------------------------------------------------

def tot_lq(xi1,xi2,lpar,lpars):

    al1,al2 = call_lq_nie(xi1,xi2,lpar)
    for i in lpars:
        al1s,al2s = call_lq_nie(xi1,xi2,i)
        al1 = al1+al1s
        al2 = al2+al2s

    yi1 = xi1-al1
    yi2 = xi2-al2

    return yi1,yi2
#@profile
#def refine_critical(lpar,lpars,critical,xi1,xi2,dsx,nfiner=4):
#    dsf = dsx/nfiner/2
#    x1t = []
#    x2t = []
#    for i in xrange(nfiner):
#        for j in xrange(nfiner):
#            x1tmp = xi1[critical>0]+(dsf*(1-nfiner)*0.5)+dsf*i
#            x2tmp = xi2[critical>0]+(dsf*(1-nfiner)*0.5)+dsf*j
#            x1t.append(x1tmp)
#            x2t.append(x2tmp)
#
#    yift1,yift2 = tot_lq(x1t,x2t,lpar,lpars)
#
#    return yift1,yift2
def refine_critical(lpar,lpars,critical,xi1,xi2,dsx,nfiner=4):
    dsf = dsx/nfiner/4
    x1t = []
    x2t = []
    x1 = xi1[critical>0]
    x2 = xi2[critical>0]
    for i in xrange(nfiner):
        for j in xrange(nfiner):
            x1tmp = x1+(dsf*(1-nfiner)*0.5)+dsf*i
            x2tmp = x2+(dsf*(1-nfiner)*0.5)+dsf*j
            x1t.append(x1tmp)
            x2t.append(x2tmp)

    yift1,yift2 = tot_lq(x1t,x2t,lpar,lpars)

    return yift1,yift2

@profile
def lensed_images(xi1,xi2,spar,lpar,lpars):

    dsx = xi1[1,1]-xi1[0,0]
    al1,al2 = call_lq_nie(xi1,xi2,lpar)
    for i in lpars:
        al1s,al2s = call_lq_nie(xi1,xi2,i)
        al1 = al1+al1s
        al2 = al2+al2s

    a12,a11 = np.gradient(al1,dsx)
    a22,a21 = np.gradient(al2,dsx)

    mu = 1.0/(1.0-(a11+a22)+a11*a22-a12*a21)

    s_image = call_gauss_2d(xi1,xi2,spar)

    yi1 = xi1-al1
    yi2 = xi2-al2

    g_lensimage = call_gauss_2d(yi1,yi2,spar)

    return s_image,g_lensimage,mu,yi1,yi2

#--------------------------------------------------------------------
def lens_images(xi1,xi2,gpar,gpars):

    g_lens = call_gauss_2d(xi1,xi2,gpar)
    for i in gpars:
        g_lens_subs = call_gauss_2d(xi1,xi2,i)
        g_lens = g_lens + g_lens_subs
    return g_lens

def mmbr_images(xi1,xi2,gpar,gpars):

    g_lens = call_tophat_2d(xi1,xi2,gpar)
    g_edge = call_find_critical_curve(g_lens)
    for i in gpars:
        g_lens_subs = call_tophat_2d(xi1,xi2,i)
        g_edge_subs = call_find_critical_curve(g_lens_subs)
        g_edge = g_edge+g_edge_subs
    g_edge[g_edge>0.0] = 1.0
    return g_edge

def main():
    nnn = 1024
    boxsize = 4.0
    dsx = boxsize/nnn
    xi1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi1,xi2 = np.meshgrid(xi1,xi2)


    baset = np.zeros((nnn,nnn,3),'uint8')
    base0 = np.zeros((nnn,nnn,3),'uint8')
    base1 = np.zeros((nnn,nnn,3),'uint8')
    base2 = np.zeros((nnn,nnn,3),'uint8')
    base3 = np.zeros((nnn,nnn,3),'uint8')
    base4 = np.zeros((nnn,nnn,3),'uint8')

    #----------------------------------------------------
    # parameters of source
    x = 0
    y = 0
    #step = 1
    gr_sig = 0.02
    gr_eq = 1.0
    gr_pa = 0.0

    #----------------------------------------------------
    # lens parameters for mainhalo
    xlc0 = 0.0
    ylc0 = 0.0
    ql0 = 0.7
    rc0 = 0.1
    re0 = 1.0
    phi0 = 0.0
    #----------------------------------------------------
    # lens parameters for subhalo
    xlcs = 0.7
    ylcs = 0.77
    qls = 0.999999
    rcs = 0.000001
    res = 0.05
    phis = 0.0
    #----------------------------------------------------
    # parameters of NIE model (lens model, deflection angles)
    #----------------------------------------------------
    # 1, y position of center
    # 2, x position of center
    # 3, minor-to-major axis ratio
    # 4, size of flat core
    # 5, Einstein radius (lensing strength)
    # 6, major-axis position angle (degrees) c.c.w. from y axis
    lpar_sub = np.asarray([ylcs,xlcs,qls,rcs,res,phis])
    lpars = [lpar_sub]

    #----------------------------------------------------
    # luminosity parameters for mainbhalo
    ap0 = 1.0
    l_sig0 = 0.5
    #----------------------------------------------------
    # luminosity parameters for subhalo
    aps = 0.4
    l_sigs = 0.05
    #----------------------------------------------------
    # Parameters of Gaussian model (luminosity distribution of lenses)
    #----------------------------------------------------
    # 1, peak brightness value
    # 2, Gaussian "sigma" (i.e., size)
    # 3, y position of center
    # 4, x position of center
    # 5, minor-to-major axis ratio
    # 6, major-axis position angle (degrees) c.c.w. from y axis

    gpars_sub = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis])
    gpars = [gpars_sub]
    #---------------------------------------------------

    lpar =  np.asarray([ylc0,xlc0,ql0,rc0,re0,phi0])
    gpar = np.asarray([ylc0,xlc0,ql0,ap0,l_sig0,phi0])

    #----------------------------------------------
    #parameters of source galaxies.
    #----------------------------------------------
    g_amp = 1.0         # peak brightness value
    g_sig = gr_sig          # Gaussian "sigma" (i.e., size)
    g_ycen = y*2.0/nnn  # y position of center
    g_xcen = x*2.0/nnn  # x position of center
    g_axrat = gr_eq       # minor-to-major axis ratio
    g_pa = gr_pa          # major-axis position angle (degrees) c.c.w. from y axis
    spar = np.asarray([g_ycen,g_xcen,g_axrat,g_amp,g_sig,g_pa])
    #----------------------------------------------

    g_lenses = lens_images(xi1,xi2,gpar,gpars)
    g_shapes = mmbr_images(xi1,xi2,gpar,gpars)

    baset[:,:,0] = g_shapes*255
    baset[:,:,1] = g_shapes*255
    baset[:,:,2] = g_shapes*255

    s_image,g_lensimage,mu,yi1,yi2 = lensed_images(xi1,xi2,spar,lpar,lpars)
    mu = 1.0/mu

    base0[:,:,0] = g_lenses*255
    base0[:,:,1] = g_lenses*127
    base0[:,:,2] = g_lenses*0

    base1[:,:,0] = s_image*255
    base1[:,:,1] = s_image*255
    base1[:,:,2] = s_image*255

    base2[:,:,0] = g_lensimage*102
    base2[:,:,1] = g_lensimage*178
    base2[:,:,2] = g_lensimage*255

    critical = call_find_critical_curve(mu)
    base3[:,:,0] = critical*255
    base3[:,:,1] = critical*0
    base3[:,:,2] = critical*0

    yif1,yif2 = refine_critical(lpar,lpars,critical,xi1,xi2,dsx)
    caustic = call_forward_cic(nnn,nnn,boxsize,yif1.flat,yif2.flat)
    caustic[caustic>0]=1

    base4[:,:,0] = caustic*0
    base4[:,:,1] = caustic*255
    base4[:,:,2] = caustic*0

    wf = base1+base2+base3+base4
    print np.max(wf)
    #pl.figure()
    #pl.contourf(wf[:,:,0])
    #pl.figure()
    #pl.contourf(wf[:,:,1])
    #pl.figure()
    #pl.contourf(wf[:,:,2])
    #pl.show()

if __name__ == '__main__':
    main()
