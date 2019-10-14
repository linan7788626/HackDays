#!/usr/bin/env python
import numpy as np
import libv4_cv as lv4
import pyfits
import subprocess as sp
import mycosmology as mm
import sys

def re0_sigma(sigma,zl,zs):
    res = 4.0*np.pi*(sigma/mm.cv)**2.0*mm.Da2(zl,zs)/mm.Da2(zs)
    return res

def hfunc(x1,x2,rcore,qe):
    res = np.sqrt(qe*qe*(rcore*rcore+x1*x1)+x2*x2)
    return res


def nie_phi(x1,x2,re0,rcore,qe):

    res0 = re0/np.sqrt(1-qe*qe)
    al1 = res0*np.arctan(x1*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore))
    al2 = res0*np.arctanh(x2*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore*qe*qe))

    res1 = x1*al1+x2*al2
    res2 = re0*rcore*np.log(np.sqrt((hfunc(x1,x2)+rcore)**2.0+(1-qe*qe)*x1*x1))
    res = res1-res2
    return res

def nie_alphas(xi1,xi2,xc1,xc2,re0,rcore,qe):
    x1 = xi1-xc1
    x2 = xi2-xc1
    res0 = re0/np.sqrt(1-qe*qe)
    al1 = np.arctan(x1*np.sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore))
    al2 = np.arctanh(x2*np.sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore*qe*qe))
    return res0*al1,res0*al2

def nie_kappa(x1,x2,re0,rcore,qe):
    res = re0/(2.0*np.sqrt(qe*qe*(rcore+x1*x1)+x2*x2))
    return res

def nie_mu(x1,x2,re0,rcore,qe):
    res = 1.0/(1.0-re0/hfunc(x1,x2)-re0*re0*rcore/(hfunc(x1,x2)*((hfunc(x1,x2)+rcore)**2+(1-qe*qe)*x1*x1)))
    return res

def new_nie_all(xi1,xi2,lpar):
    xc1 = lpar[0]
    xc2 = lpar[1]
    b = lpar[2]
    s = lpar[3]
    q = lpar[4]
    rot = lpar[5]

    dsx = xi1[1,1]-xi1[0,0]

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    a1 = b/np.sqrt(1-q*q)*np.arctan(x1*np.sqrt(1-q*q)/(wx+s))
    a2 = b/np.sqrt(1-q*q)*np.arctanh(x2*np.sqrt(1-q*q)/(wx+q*q*s))

    hx = np.sqrt((wx+s)**2.0+(1-q*q)*x1*x1)
    phi = x1*a1+x2*a2-b*s*np.log(hx)+b*q*s*np.log((1+q)*s)

    ai2,ai1 = np.gradient(phi,dsx)

    return phi,ai1,ai2#,kappa,mu,y1,y2

def nie_all(xi1,xi2,xc1,xc2,b,s,q,rot,ys1,ys2):

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    al1 = b/np.sqrt(1-q*q)*np.arctan(x1*np.sqrt(1-q*q)/(wx+s))
    al2 = b/np.sqrt(1-q*q)*np.arctanh(x2*np.sqrt(1-q*q)/(wx+q*q*s))

    kappa = b/(2.0*wx)

    hx = np.sqrt((wx+s)**2.0+(1-q*q)*x1*x1)
    phi = x1*al1+x2*al2-b*s*np.log(hx)+b*q*s*np.log((1+q)*s)

    Kc = 1.0
    #Kc = (1.0+zl)/c*(Dl*Ds/Dls)
    td = Kc*(0.5*((al1)**2.0+(al2)**2.0)-phi)
    #td = Kc*(0.5*((x1-ys1)**2.0+(x2-ys2)**2.0)-phi)

    y1 = x1-al1
    y2 = x2-al2

    y1,y2 = xy_rotate(y1,y2,xc1,xc2,-rot)

#------------------------------------------------------------------
    demon1 = ((wx+s)**2+(1.0-q*q)*x1*x1)*wx
    demon2 = (((wx+q*q*s)**2-(1.0-q*q)*x2*x2)*wx)
    y11 = 1-b*(wx*(wx+s)-q*q*x1*x1)/demon1
    y22 = 1-b*(wx*(wx+q*q*s)-x2*x2)/demon2
    y12 = -b*x1*x2/demon1
    y21 = -b*x1*x2*q*q/demon2

    mu = 1.0/(y11*y22-y12*y21)

    return phi,td,al1,al2,kappa,mu,y1,y2

def multiple_nie_all(xi1,xi2,lpars_list):
    phi = xi1*0.0
    al1 = xi1*0.0
    al2 = xi1*0.0
    for i in lpars_list:
        phi_tmp,al1_tmp,al2_tmp = lpar_nie_all(xi1,xi2,i)
        phi = phi + phi_tmp
        al1 = al1 + al1_tmp
        al2 = al2 + al2_tmp

    return phi,al1,al2

def multiple_new_nie_all(xi1,xi2,lpars_list):
    phi = xi1*0.0
    al1 = xi1*0.0
    al2 = xi1*0.0
    for i in lpars_list:
        phi_tmp,al1_tmp,al2_tmp = new_nie_all(xi1,xi2,i)
        phi = phi + phi_tmp
        al1 = al1 + al1_tmp
        al2 = al2 + al2_tmp

    return phi,al1,al2

def lpar_nie_all(xi1,xi2,lpar):

    xc1 = lpar[0]
    xc2 = lpar[1]
    b = lpar[2]
    s = lpar[3]
    q = lpar[4]
    rot = lpar[5]

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    al1 = b/np.sqrt(1-q*q)*np.arctan(x1*np.sqrt(1-q*q)/(wx+s))
    al2 = b/np.sqrt(1-q*q)*np.arctanh(x2*np.sqrt(1-q*q)/(wx+q*q*s))

    hx = np.sqrt((wx+s)**2.0+(1-q*q)*x1*x1)
    phi = x1*al1+x2*al2-b*s*np.log(hx)+b*q*s*np.log((1+q)*s)

    return phi,al1,al2

def lensed_images(xi1,xi2,yi1,yi2,gpar):

    g_image = gauss_2d(xi1,xi2,gpar)
    g_lensimage = gauss_2d(yi1,yi2,gpar)

    return g_image,g_lensimage

def lensed_images_point(xi1,xi2,yi1,yi2,gpar):

    g_image = tophat_2d(xi1,xi2,gpar)
    g_lensimage = tophat_2d(yi1,yi2,gpar)

    return g_image,g_lensimage

def xy_rotate(x, y, xcen, ycen, phi):

    phirad = np.deg2rad(phi)
    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
    return (xnew,ynew)

def tophat_2d(x, y, par):

    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4])/np.abs(par[1])**2
    res = r_ell_sq*0.0
    res[r_ell_sq<=1.0] = par[0]
    return res

def gauss_2d(x, y, par):

    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / np.abs(par[1])**2
    return par[0] * np.exp(-0.5*r_ell_sq)

def gauss_1d(x, x0,sigma,a):

    r_ell_sq = ((x-x0)**2.0/sigma**2)
    res = a*np.exp(-0.5*r_ell_sq)
    return res

def parabola_1d(x,xb,xc,a):
    res = a*(x-xb)*(2.0*xc-(x-xb))
    res[res<=0] = 0.0
    return res

def find_critical_curve(mu):
    rows,cols = np.indices(np.shape(mu))
    cdtn = np.sign(mu)*(np.sign(mu[rows-1,cols])+np.sign(mu[rows,cols-1])+np.sign(mu[(rows+1)%len(rows),cols])+np.sign(mu[rows,(cols+1)%len(cols)]))

    res = mu*0
    res[cdtn<4] = 1
    res[cdtn>=4] = 0

    return res

def lens_galaxies(xi1,xi2,glpar):

    g_lens = gauss_2d(xi1,xi2,glpar)

    return g_lens

def main():

    nnn = 512
    dsi = 0.2
    dsx = 0.2
    boxsize = 0.2*nnn
    zzl = 0.5
    zs0 = 10.0

    xi1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi1,xi2 = np.meshgrid(xi1,xi2)

    limages = np.zeros((nnn,nnn))
    #----------------------------------------------------
    # lens parameters for main halo
    xlc1 = -1.5
    xlc2 = 0.0
    ql0 = 0.699999999999
    rc0 = 0.000000000001
    re0 = 6.0
    #phi0 = 30.0
    band = sys.argv[1]

    ai01,ai02 = nie_alphas(xi1,xi2,xlc1,xlc2,re0,rc0,ql0)
    gimage = pyfits.getdata("../gals_lens_galsim/test_"+band+".fits")
    #----------------------------------------------------


    ysc1 = 0.0
    ysc2 = 0.0
    factor_z = 0;
    factor_z0 = mm.Da(zs0)/mm.Da2(zzl,zs0);

    s = " "
    input_fits_dir = "../catsim_galsim_fits/"
    input_fits_sh = s.join(("ls", input_fits_dir,"|grep",band+".fits"))
    input_fits = sp.check_output(input_fits_sh,shell=True)
    input_fits = input_fits.rstrip()
    input_fits_list = input_fits.split("\n")

    for i in input_fits_list:

        array = i.split("_")
        zzsd = np.double(array[1])
        zzsu = np.double(array[2])
        zzs = (zzsd+zzsu)*0.5

        srcs = pyfits.getdata(input_fits_dir+i)
        srcs = np.array(srcs,dtype="<d")

        if zzs < zzl:
            factor_z=0.0
        else:
            factor_z=mm.Da2(zzl,zzs)/mm.Da(zzs)*factor_z0

        #print factor_z

        ai1 = ai01*factor_z
        ai2 = ai02*factor_z

        yi1 = xi1 - ai1
        yi2 = xi2 - ai2

        #g_sn_pin = lv4.call_ray_tracing(g_sn,xi1,xi2,ysc1,ysc2,dsi)
        limages_tmp = lv4.call_ray_tracing(srcs,yi1,yi2,ysc1,ysc2,dsi)
        limages = limages + limages_tmp

    #limages[128:384,128:384] = limages[128:384,128:384]+gimage
    output_file_name = "".join(("../lensed_images/lensed_twinkles_",band,".fits"))
    pyfits.writeto(output_file_name,limages,clobber=True)

    limages[128:384,128:384] = limages[128:384,128:384]+gimage
    output_file_name = "".join(("../final_images/final_twinkles_",band,".fits"))
    pyfits.writeto(output_file_name,limages,clobber=True)

    return 0


if __name__ == '__main__':
    main()

