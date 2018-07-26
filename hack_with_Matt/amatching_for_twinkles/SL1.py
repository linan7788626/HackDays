#!/usr/bin/env python
import numpy as np
import os.path
from astropy.cosmology import Planck13 as p13
from astropy import constants as const
import scipy.special as spf
import pylab as pl
import astropy.io.fits as pyfits

apr = 206269.43 # arcseconds per rad
name="lens.cat"
m = np.loadtxt(name, skiprows=0)

    #--------------------------------------------------------------------
    # Construct regular grids
    #
def make_r_coor(nc,dsx):

    bsz = nc*dsx
    x1 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0

    x2,x1 = np.meshgrid(x1,x2)
    return x1,x2
#--------------------------------------------------------------------
# Calculate deflection angles
#
def lens_equation_sie(x1,x2,xc1,xc2,q,rc,re,pha):

    phirad = np.deg2rad(pha)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1-xc1)*cosa+(x2-xc2)*sina
    xt2 = (x2-xc2)*cosa-(x1-xc1)*sina

    phi = np.sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc)
    sq = np.sqrt(1.0-q*q)
    pd1 = phi+rc/q
    pd2 = phi+rc*q
    fx1 = sq*xt1/pd1
    fx2 = sq*xt2/pd2
    qs = np.sqrt(q)

    a1 = qs/sq*np.arctan(fx1)
    a2 = qs/sq*np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 =-sina

    fx11 = xt11/pd1-xt1*(xt1*q*q*xt11+xt2*xt21)/(phi*pd1*pd1)
    fx22 = xt22/pd2-xt2*(xt1*q*q*xt12+xt2*xt22)/(phi*pd2*pd2)
    fx12 = xt12/pd1-xt1*(xt1*q*q*xt12+xt2*xt22)/(phi*pd1*pd1)
    fx21 = xt21/pd2-xt2*(xt1*q*q*xt11+xt2*xt21)/(phi*pd2*pd2)

    a11 = qs/(1.0+fx1*fx1)*fx11
    a22 = qs/(1.0-fx2*fx2)*fx22
    a12 = qs/(1.0+fx1*fx1)*fx12
    a21 = qs/(1.0-fx2*fx2)*fx21

    rea11 = (a11*cosa-a21*sina)*re
    rea22 = (a22*cosa+a12*sina)*re
    rea12 = (a12*cosa-a22*sina)*re
    rea21 = (a21*cosa+a11*sina)*re

    y11 = 1.0-rea11
    y22 = 1.0-rea22
    y12 = 0.0-rea12
    y21 = 0.0-rea21

    jacobian = y11*y22-y12*y21
    mu = 1.0/jacobian

    res1 = (a1*cosa-a2*sina)*re
    res2 = (a2*cosa+a1*sina)*re
    return res1,res2,mu
#--------------------------------------------------------------------
# Rotate regular grids
#
def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x-xcen)*np.cos(phirad)+(y-ycen)*np.sin(phirad)
    ynew = (y-ycen)*np.cos(phirad)-(x-xcen)*np.sin(phirad)
    return (xnew,ynew)
#--------------------------------------------------------------------
# Calculate Einstein Radius according to Velocity Dispersion
#
def re_sv(sv,z1,z2):
    Da_s = p13.angular_diameter_distance(z2).value
    Da_ls = p13.angular_diameter_distance_z1z2(z1,z2).value
    res = 4.0*np.pi*(sv**2.0/(const.c.value/1e3)**2.0)*Da_ls/Da_s*apr
    return res
#--------------------------------------------------------------------
# 2D Sersic Profile
#
def sersic_2d(xi1,xi2,xc1,xc2,Ieff,Reff_arc,ql,pha,ndex=4.0):
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_scale = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff_arc
    res = Ieff*np.exp(-bn*((R_scale)**(1.0/ndex)-1.0))
    return res
#--------------------------------------------------------------------
# Convert total magnitude to effective flux
#
def sersic_mag_tot_to_Ieff(mag_tot,Reff,ndex,z,Rcut=100): # g band
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    xtmp = bn*(Rcut)**(1.0/ndex)
    ktmp = 2.0*np.pi*ndex*(np.e**bn/(bn**(2.0*ndex)))*spf.gammainc(2.0*ndex,xtmp)*spf.gamma(2.0*ndex)
    Dl_s = p13.luminosity_distance(z).value*1e6

    Ieff = 10.0**((5.12-2.5*np.log10(ktmp)-5.0*np.log10(Reff)+5.0*np.log10(Dl_s/10)-mag_tot)*0.4) # L_sun/pc^2
    return Ieff
#--------------------------------------------------------------------
# Generate lensed Sersic Profile
#
def lensed_sersic(xi1,xi2,source_cat,lens_cat):

    xlc1 = lens_cat[0]          # x position of the lens, arcseconds
    xlc2 = lens_cat[1]          # y position of the lens, arcseconds
    rlc  = 0.0                  # core size of Non-singular Isothermal Ellipsoid
    rle  = re_sv(lens_cat[2],lens_cat[5],source_cat[7])       #Einstein radius of lens, arcseconds.
    ql   = lens_cat[3]          # axis ratio b/a
    phl = lens_cat[4]           # orientation, degree
    #----------------------------------------------------------------------
    ai1,ai2,mua = lens_equation_sie(xi1,xi2,xlc1,xlc2,ql,rlc,rle,phl)

    yi1 = xi1-ai1
    yi2 = xi2-ai2

    ysc1 = source_cat[0]        # x position of the source, arcseconds
    ysc2 = source_cat[1]        # y position of the source, arcseconds
    mag_tot = source_cat[2]     # total magnitude of the source
    Reff_arc = source_cat[3]    # Effective Radius of the source, arcseconds
    qs   = source_cat[4]        # axis ratio of the source, b/a
    phs  = source_cat[5]        # orientation of the source, degree
    ndex = source_cat[6]        # index of the source
    zs = source_cat[7]          # redshift of the source

    Da_s = p13.angular_diameter_distance(zs).value
    Reff = Reff_arc*Da_s/apr*1e6 # pc
    Ieff = sersic_mag_tot_to_Ieff(mag_tot,Reff,ndex,zs) #L_sun/kpc^2
    g_limage = sersic_2d(yi1,yi2,ysc1,ysc2,Ieff,Reff_arc,qs,phs,ndex)
    g_source = sersic_2d(xi1,xi2,ysc1,ysc2,Ieff,Reff_arc,qs,phs,ndex)

    mag_lensed = mag_tot - np.log(np.sum(g_limage)/np.sum(g_source))

    return mag_lensed, g_limage

if __name__ == '__main__':
    for p in range(0,2):
        print p
    
    #-----------------------------------------------------------------------
    # parameters of the mass model of the lens
    #
        xl1 = m[p,0]
        xl2 = m[p,1]
        vd =  m[p,2]      
        ql  = m[p,3]     
        phl = m[p,4]       
        zl = m[p,5] 

        lens_cats = [xl1, xl2, vd, ql, phl, zl]
        print lens_cats
        #-----------------------------------------------------------------------
        # parameters of the host galaxy of lensed point source
        #

        ysc1 = m[p,6]      # projected postion of the source, arcsecond
        ysc2 = m[p,7]      # projected postion of the source, arcsecnd
        mags = m[p,8]      # magnitude of the source
        Reff = m[p,9]      # Effective radius of Sersic Profile of the source, arcsecnd
        qs = m[p,10]        # b/a of the source
        phs = m[p,11]      # orientation of the source, degree
        ns = m[p,12]          # index of the Sersic Profile of the source
        zs = m[p,13]         # redshift of the source

        srcs_cats = [ysc1, ysc2, mags, Reff, qs, phs, ns, zs]
        print srcs_cats
        #-----------------------------------------------------------------------
        # Calculate the lensed images
        #
        dsx = 0.01  # pixel size per side, arcseconds
        nnn = 1000  # number of pixels per side
        xi1, xi2 = make_r_coor(nnn,dsx)
        #-----------------------------------------------------------------------
        # Calculate the lensed image and its total magnitude
        #
        lensed_mag, lensed_image = lensed_sersic(xi1,xi2,srcs_cats,lens_cats)

        f = open('workfile', 'w')

        if p == 0:
            phile = open("SL.out","w")
        else:
            phile = open("SL.out","a")

        print >> phile, lensed_mag

        if os.path.isfile('lenser'+str(p)+'.fits'):    
            os.remove('lenser'+str(p)+'.fits')
            hdu = pyfits.PrimaryHDU(lensed_image)
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('lenser'+str(p)+'.fits')
        else:  
            hdu = pyfits.PrimaryHDU(lensed_image)
            hdulist = pyfits.HDUList([hdu])
            hdulist.writeto('lenser'+str(p)+'.fits')          
            
        phile = open('lenser'+str(p)+'.out','w')
        print >>phile, lensed_image[0], lensed_image[1], lensed_image[2], lensed_image[3], lensed_image[4], lensed_image[5], lensed_image[6], lensed_image[7]
       
        phile.close()

    #-----------------------------------------------------------------------
    # Visualize the lensed imagesr. If we normalize the image, the 2D matrix
    # becomes the normalized histogram of the light distribution.
    #
    #pl.contourf(lensed_image)
    #pl.colorbar()
    #pl.show()
