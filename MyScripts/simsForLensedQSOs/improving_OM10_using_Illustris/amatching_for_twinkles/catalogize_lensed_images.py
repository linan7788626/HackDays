#!/usr/bin/env python
import numpy as np
from astropy.cosmology import Planck13 as p13
from astropy import constants as const
import pylab as pl
import scipy.special as spf
import om10
import pyfits

apr = 206269.43
# m = M + 5 log10 d/(1pc) - 5
#--------------------------------------------------------------------
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

#--------------------------------------------------------------------
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
def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x-xcen)*np.cos(phirad)+(y-ycen)*np.sin(phirad)
    ynew = (y-ycen)*np.cos(phirad)-(x-xcen)*np.sin(phirad)
    return (xnew,ynew)

def gauss_2d(x, y, par):
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    R_arc = np.sqrt(((xnew**2)*par[4]+(ynew**2)/par[4]))/np.abs(par[1])
    res = par[0]*np.exp(-R_arc**2.0)
    return res

def re_sv(sv,z1,z2):
    Da_s = p13.angular_diameter_distance(z2).value
    Da_ls = p13.angular_diameter_distance_z1z2(z1,z2).value
    res = 4.0*np.pi*(sv**2.0/(const.c.value/1e3)**2.0)*Da_ls/Da_s*apr
    return res

def I2HSTccd(input_I,dsx):
    #dl = mm.Dl(zl)
    #dsp = dsx*mm.Dl(zl)/mm.apr

    #mag = 5.12-2.5*np.log(input_I*dsp*dsp)+5.0*np.log10(dl*1e6/10.0) # mag per arcsecnd^2
    #mag = 5.12+5.0*np.log10(mm.apr/dsx)-2.5*np.log10(input_I*100.0)
    mag = 5.12+(5.0*np.log10(mm.apr/dsx)-5)-2.5*np.log10(input_I)

    # Zero point of ACS F606w in AB mag
    zp_candels_ABmag  =26.49113

    # Flux of a given magnitude (input_mag), erg cm-2 s-1 Ang-1
    Flux = 10.0**(-0.4*(mag-zp_candels_ABmag))

    # PHOTFLAM: inverse sensitivity, erg cm-2 s-1 Ang-1
    #PhotFlam = 7.862481e-20

    # Pixel value, Number of count per second
    Ncount = Flux#/PhotFlam

    return Ncount

def Brightness(Re,Vd):
    A = 2.09
    B = -1.49
    C = - 1.19

    log10_Ie_bar   =A*np.log10(Vd)+B*np.log10(Re)+C # Kpc, km/s
    I_e = 10**(log10_Ie_bar)/3.61 # http://www.astro.spbu.ru/staff/resh/Lectures/lec2.pdf
    #mag_e_bar = 5.14+21.572-log_Ie_bar/0.4 # mag per arcsecnd^2
    #mag_e = mag_e_bar+1.39

    return I_e

def de_vaucouleurs_2d(xi1, xi2, xc1, xc2, I_e, R_arc, ql, pha):
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_arc = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/R_arc
    #res = mu_e + 8.32678*((R_arc)**0.25-1.0)
    res = I_e*10.0**(-3.33071*((R_arc)**0.25-1.0))
    res_th = I_e*10.0**(-3.33071*((0.03)**0.25-1.0))
    res[res>res_th] = res_th

    return res

def sersic_classic_2d(xi1,xi2,xc1,xc2,Ieff,Reff,ql,pha,ndex=4.0):
    #bn = 1.9992*ndex-0.3271 #(0.5<ndex<10.0)
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_scale = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff
    res = Ieff*np.exp(-bn*((R_scale)**(1.0/ndex)-1.0))
    return res

def sersic_2d(xi1,xi2,xc1,xc2,Ieff,Reff_arc,ql,pha,ndex=4.0):
    #bn = 1.9992*ndex-0.3271 #(0.5<ndex<10.0)
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_scale = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff_arc
    res = Ieff*np.exp(-bn*((R_scale)**(1.0/ndex)-1.0))
    return res

def sersic_L_Rcut(Ieff,Reff,ndex,Rcut=100):
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    xtmp = bn*(Rcut)**(1.0/ndex)
    res = 2.0*np.pi*Ieff*Reff**2.0*ndex*(np.e**bn/(bn**(2.0*ndex)))*spf.gammainc(2.0*ndex,xtmp)*spf.gamma(2.0*ndex)

    #res = np.pi*Ieff*Reff**2.0*(2.0*ndex)*(np.e**bn/(bn**(2.0*ndex)))*spf.gamma(2.0*ndex)
    return res

def sersic_Iebar_to_Ieff(Iebar,Reff,ndex):
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    res = Iebar/(2.0*ndex*(np.e**bn/(bn**(2.0*ndex)))*spf.gammainc(2.0*ndex,bn)*spf.gamma(2.0*ndex))
    print (2.0*ndex*(np.e**bn/(bn**(2.0*ndex)))*spf.gammainc(2.0*ndex,bn)*spf.gamma(2.0*ndex))
    return res

def sersic_Abs_MAG_Rcut(Ieff,Reff,ndex):
    LB = sersic_L_Rcut(Ieff,Reff,ndex)
    res = 5.12 - 2.5*np.log10(LB)
    return res


def sersic_mag_in_Rcut(Ieff,Reff,ndex,z,Rcut=100): # g band
    LB = sersic_L_Rcut(Ieff,Reff,ndex)
    Dl_l = p13.luminosity_distance(z).value*1e6
    res = 5.12 - 2.5*np.log10(LB) + 5.0*np.log10(Dl_l/10)
    return res

def sersic_mag_tot_to_Ieff(mag_tot,Reff,ndex,z,Rcut=100): # g band
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    xtmp = bn*(Rcut)**(1.0/ndex)
    ktmp = 2.0*np.pi*ndex*(np.e**bn/(bn**(2.0*ndex)))*spf.gammainc(2.0*ndex,xtmp)*spf.gamma(2.0*ndex)
    Dl_s = p13.luminosity_distance(z).value*1e6

    Ieff = 10.0**((5.12-2.5*np.log10(ktmp)-5.0*np.log10(Reff)+5.0*np.log10(Dl_s/10)-mag_tot)*0.4) # L_sun/pc^2
    return Ieff

def sersic_mag_tot_to_mueff(mag_tot,Reff,ndex,z,Rcut=100): # g band
    Ieff = sersic_mag_tot_to_Ieff(mag_tot,Reff,ndex,z)
    mueff = 5.12+21.572-2.5*np.log10(Ieff)
    return mueff

def sersic_SB_2D_tot(xi1,xi2,xc1,xc2,mag_tot,Reff_arc,ql,pha,z,ndex=4.0):

    #Dl_s = p13.luminosity_distance(z).value
    Da_s = p13.angular_diameter_distance(z).value
    Reff = Reff_arc*Da_s/apr*1e6
    mueff=sersic_mag_tot_to_mueff(mag_tot,Reff,ndex,z)
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_arc = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff_arc
    res = mueff + (2.5*bn)/np.log(10.0)*((R_arc/Reff_arc)**(1.0/ndex)-1.0)
    return res

def sersic_SB_2D(xi1,xi2,xc1,xc2,mueff,Reff_arc,ql,pha,ndex=4.0):
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_arc = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff_arc
    res = mueff + (2.5*bn)/np.log(10.0)*((R_arc/Reff_arc)**(1.0/ndex)-1.0)
    return res

def mag_tot_SED(band_in,mag_tot_in, band_out, mag_tot_out):
    #####
    ####
    ####
    #####
    return mag_tot_out

def sersic_I0(Stotal,Reff,ndex): # at Reff
    res = Stotal(2.0*np.pi*Reff**2.0*ndex*spf.gamma(2.0*ndex))
    return res

def producing_lensed_images(xi1,xi2,source_cat,lens_cat):

    xlc1 = lens_cat[0]
    xlc2 = lens_cat[1]
    ql   = lens_cat[2]
    rlc  = 0.0
    rle  = re_sv(lens_cat[3],lens_cat[5],source_cat[7])       #Einstein radius of lens Arc Seconds.

    phl = lens_cat[4]
    #----------------------------------------------------------------------
    ai1,ai2,mua = lens_equation_sie(xi1,xi2,xlc1,xlc2,ql,rlc,rle,phl)

    yi1 = xi1-ai1*0.0
    yi2 = xi2-ai2*0.0

    ysc1 = source_cat[0]
    ysc2 = source_cat[1]
    mag_tot = source_cat[2]
    Reff_arc = np.sqrt(source_cat[3]*source_cat[4])
    qs   = np.sqrt(source_cat[3]/source_cat[4])
    phs  = source_cat[5]
    ndex = source_cat[6]
    zs = source_cat[7]

    #g_limage = sersic_SB_2D_tot(yi1,yi2,ysc1,ysc2,mag_tot,Reff_arc,qs,phs,zs)
    Da_s = p13.angular_diameter_distance(2.0).value
    Reff = Reff_arc*Da_s/apr*1e6 # pc
    Ieff = sersic_mag_tot_to_Ieff(mag_tot,Reff,ndex,zs)
    g_limage = sersic_2d(yi1,yi2,ysc1,ysc2,Ieff,Reff_arc,qs,phs,ndex)

    return g_limage

if __name__ == '__main__':

    db = om10.DB(catalog="../qso_mock.fits")
    #lid = 7176527
    lid = 630751
    lens = db.get_lens(lid)

    lid = lens.LENSID[0]
    xl1 = 0.0
    xl2 = 0.0
    vd = lens.VELDISP[0]    # needed from OM10
    zd = lens.ZLENS[0]
    zs = lens.ZSRC[0]
    ql  = 1.0 - lens.ELLIP[0]
    phi= lens.PHIE[0]
    print ql
    #how to calculate I_eff_mag, R_eff_arc, qs, phs

    #source_cat = [ysc1, ysc2, mag_tot, R_eff_arc_a, R_eff_arc_a, phs, zs]
    #lens_cat = [xc1, xc2, ql, sigma_v, phl, zl]
    mag_srcs = 25

    ndex = 4

    ys1 = 0.1
    ys2 = 0.2
    resa = 0.4
    resb = 0.2
    stheta = 35.0

    source_cat = [ys1, ys2, mag_srcs, resa, resb, stheta, ndex, zs]
    lens_cat = [xl1, xl2, ql, vd, phi, zd, lid]

    dsx = 0.01
    nnn = 1000

    xi1, xi2 = make_r_coor(nnn,dsx)

    lensed_images = producing_lensed_images(xi1,xi2,source_cat,lens_cat)


    Dl_s = p13.luminosity_distance(zs).value
    Da_s = p13.angular_diameter_distance(zs).value
    Reff_arc = np.sqrt(resa*resb)
    Reff = Reff_arc*Da_s/apr*1e6 # pc
    ndex = 4.0
    Ieff = sersic_mag_tot_to_Ieff(mag_srcs,Reff,ndex,zs)
    print Reff, Ieff
    Ltot = sersic_L_Rcut(Ieff,Reff,ndex)

    MAG = sersic_Abs_MAG_Rcut(Ieff,Reff,ndex)
    print MAG
    mueff = sersic_mag_tot_to_mueff(24,Reff,ndex,2.0)
    Ieff2 = 10.0**((21.572+5.12-mueff)*0.4)
    #Ltmp = np.sum(sersic_2d(xi1,xi2,0.0,0.0,Ieff,Reff_arc,0.5,35.0,ndex))*(dsx*Da_s/apr*1e6)**2.0
    Ltmp = np.sum(lensed_images)*(dsx*Da_s/apr*1e6)**2.0
    magtmp = 5.12 - 2.5*np.log10(Ltmp)+5.0*np.log10(Dl_s*1e6/10.0)

    #print magtmp

    pyfits.writeto(str(magtmp)+"_lensed_images.fits",lensed_images, clobber=True)


    #pl.figure()
    #pl.contourf(xi1,xi2,lensed_images)
    #pl.colorbar()
    #pl.show()
