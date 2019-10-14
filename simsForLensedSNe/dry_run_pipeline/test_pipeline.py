#!/usr/bin/env python
import numpy as np
import libv4_cv as lv4
import mycosmology as mm
import scipy.signal as ss
import astropy.io.fits as pyfits
from astropy.cosmology import Planck13
import pylab as pl
import scipy.interpolate as sci
import pixcos2pixsdss as p2p
import congrid

def mags_to_ncounts(input_array,zeropoint):
    res = input_array
    return res

def ncounts_to_flux(input_array,zeropoint):
    res = input_array
    return res

def mag_to_flux(input_array,zeropoint):
    res = input_array
    return res

def rebin_psf(input_psf,new_shape):
    nxo,nyo = np.shape(input_psf)
    nxn,nyn = new_shape

    xo = np.linspace(0,nxo-1.0,nxo)+0.5
    yo = np.linspace(0,nyo-1.0,nyo)+0.5
    xo,yo = np.meshgrid(xo,yo)
    xo=xo.reshape((nxo*nyo))
    yo=yo.reshape((nxo*nyo))
    zo=input_psf.reshape((nxo*nyo))


    xn = np.linspace(0,nxo-1.0,nxn)+0.5
    yn = np.linspace(0,nyo-1.0,nyn)+0.5
    xn,yn = np.meshgrid(xn,yn)

    print np.max(xo),np.min(xo)
    print np.max(xn),np.min(xn)

    res = sci.griddata(np.array([xo,yo]).T, zo, (xn, yn), method='linear')
    return res


#def re0_sigma(sigma):
#    cv = 3e5
#    Dds = 1.0
#    Ds = 2.0
#    res = 4.0*np.pi*(sigma/cv)**2.0*Dds/Ds
#    return res

nMgyCount_r=0.004760406   # nanomaggies per count for SDSS detector.
sky_r      =5.98          # SDSS typical r band sky
softbias   =1000.0        # SDSS softbias
Mgy2nanoMgy=10e+9         # nanoMaggy to Maggy
skycount=sky_r/(nMgyCount_r)


def noise_map(nx1,nx2,nstd,NoiseType):
    if NoiseType=='Poisson':
    	noise=np.random.poisson(nstd,(nx1,nx2))-nstd
    if NoiseType=='Gaussian':
    	noise=nstd*np.random.normal(0.0,1.0,(nx1,nx2))
    return noise
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
def lens_equation_sie(x1,x2,lpar):
    xc1 = lpar[0]   #x coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]   #y coordinate of the center of lens (in units of Einstein radius).
    q   = lpar[2]   #Ellipticity of lens.
    rc  = lpar[3]   #Core size of lens (in units of Einstein radius).
    re  = lpar[4]   #Einstein radius of lens.
    pha = lpar[5]   #Orintation of lens.

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
    res0 = np.sqrt(((xnew**2)*par[4]+(ynew**2)/par[4]))/np.abs(par[1])
    res = par[0]*np.exp(-res0**2.0)
    return res

def re_sv(sv,z1,z2):
    res = 4.0*np.pi*(sv**2.0/mm.vc**2.0)*mm.Da2(z1,z2)/mm.Da(z2)*mm.apr
    return res

#----Fundamental Plain------------------------------
def Brightness(Re,Vd):
    a       =1.49
    b       =0.2
    c       =-8.778
    mag_e   =((np.log10(Re)-a*np.log10(Vd)-c)/b)+20.09 # Bernardi et al 2003
    nanoMgy =Mgy2nanoMgy*10.0**(-(mag_e-22.5)/2.5)
    counts  =nanoMgy/nMgyCount_r
    return counts

def de_vaucouleurs_2d(x,y,par):
    #[I0, Re, xc1,xc2,q,pha]
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = np.sqrt((xnew**2)*par[4]+(ynew**2)/par[4])/par[1]
    #res = par[0]*np.exp(-par[1]*res0**0.25)
    res = par[0]*np.exp(-7.669*(res0**0.25-1.0))
    return res

##----de Vaucouleurs profile-------------------------
#def deVaucouleurs(x,y,xc,yc,counts,R,e,phi):
    #theta   =phi*np.pi/180.

    #xx      =x-xc
    #yy      =y-yc
    #rx      =xx*np.cos(theta)+yy*np.sin(theta)
    #ry      =-xx*np.sin(theta)+yy*np.cos(theta)
    #rr      =np.sqrt(rx*rx/(1.0-e)+ry*ry*(1.0-e))
    #image   =counts*np.exp(-7.669*((rr/R)**0.25-1.0))
    #soften  =counts*np.exp(-7.669*((0.02)**0.25-1.0))
    #ix      =np.where(image>=soften)
    #image[ix]=soften


    #return image
#--------------------------------------------------------------------
def single_run_test(ind,ysc1,ysc2,q,vd,pha,zl,zs):
    zeropoint = 18

    dsx_sdss     = 0.396         # pixel size of SDSS detector.


    R  = 2.9918     #vd is velocity dispersion.
    #zl = 0.2     #zl is the redshift of the lens galaxy.
    #zs = 1.0
    #vd = 520    #Velocity Dispersion.
    nnn = 512      #Image dimension
    bsz = 30.0 # arcsecs
    dsx = bsz/nnn         # pixel size of SDSS detector.
    nstd = 59

    xx01 = np.linspace(-bsz/2.0,bsz/2.0,nnn)+0.5*dsx
    xx02 = np.linspace(-bsz/2.0,bsz/2.0,nnn)+0.5*dsx
    xi2,xi1 = np.meshgrid(xx01,xx02)
    #----------------------------------------------------------------------
    #ysc1 = 0.2
    #ysc2 = 0.5
    dsi = 0.03
    g_source = pyfits.getdata("./439.0_149.482739_1.889989_processed.fits")
    g_source = np.array(g_source,dtype="<d")
    g_source = p2p.pixcos2pixsdss(g_source)
    #----------------------------------------------------------------------
    xc1 = 0.0       #x coordinate of the center of lens (in units of Einstein radius).
    xc2 = 0.0       #y coordinate of the center of lens (in units of Einstein radius).
    #q   = 0.7       #Ellipticity of lens.
    rc  = 0.0       #Core size of lens (in units of Einstein radius).
    re  = re_sv(vd,zl,zs)       #Einstein radius of lens.
    #pha = 45.0      #Orintation of lens.
    lpar = np.asarray([xc1,xc2,q,rc,re,pha])
    #----------------------------------------------------------------------
    ai1,ai2,mua = lens_equation_sie(xi1,xi2,lpar)

    yi1 = xi1-ai1
    yi2 = xi2-ai2

    g_limage = lv4.call_ray_tracing(g_source,yi1,yi2,ysc1,ysc2,dsi)
    g_limage = mag_to_flux(g_limage,zeropoint)

    #pl.figure()
    #pl.contourf(xi1,xi2,g_limage)
    #pl.colorbar()
    #-------------------------------------------------------------
    # Need to be Caliborate the mags
    dA = Planck13.comoving_distance(zl).value*1000./(1+zl)
    Re = dA*np.sin(R*np.pi/180./3600.)
    counts  =Brightness(R,vd)
    vpar = np.asarray([counts,Re,xc1,xc2,q,pha])
    #g_lens = deVaucouleurs(xi1,xi2,xc1,xc2,counts,R,1.0-q,pha)
    g_lens = de_vaucouleurs_2d(xi1,xi2,vpar)

    g_lens = ncounts_to_flux(g_lens*1.5e-4,zeropoint)
    #-------------------------------------------------------------
    file_psf = "../PSF_and_noise/sdsspsf.fits"
    g_psf = pyfits.getdata(file_psf)-1000.0
    g_psf = g_psf/np.sum(g_psf)
    new_shape=[0,0]
    new_shape[0]=np.shape(g_psf)[0]*dsx_sdss/dsx
    new_shape[1]=np.shape(g_psf)[1]*dsx_sdss/dsx
    g_psf = rebin_psf(g_psf,new_shape)
    print(np.max(g_psf))
    g_limage = ss.fftconvolve(g_limage+g_lens,g_psf,mode="same")

    #pl.figure()
    #pl.contourf(xi1,xi2,g_limage)
    #pl.colorbar()
    #-------------------------------------------------------------
    # Need to be Caliborate the mags
    g_noise = noise_map(nnn,nnn,nstd,"Gaussian")
    g_noise = ncounts_to_flux(g_noise*1e-0+skycount,zeropoint)
    g_limage = g_limage+g_noise

    print np.shape(g_limage)
    g_limage = congrid.congrid(g_limage,[128,128])
    g_limage = g_limage-np.min(g_limage)

    pl.figure()
    #pl.contourf(xi1,xi2,g_limage)
    pl.contourf(g_limage)
    pl.colorbar()
    #-------------------------------------------------------------

    output_filename = "../output_fits/"+str(ind)+".fits"
    pyfits.writeto(output_filename,g_limage,clobber=True)

    pl.show()

    return 0

if __name__ == '__main__':
    #from mpi4py import MPI
    #import sys
    #sourcpos = 10.0 # arcsecs
    #num_imgs = int(sys.argv[1])
    num_imgs = 1
    sourcpos = 0.0

    #comm = MPI.COMM_WORLD
    #size = comm.Get_size()
    #rank = comm.Get_rank()

    #ysc1 = np.random.random(num_imgs)*sourcpos-sourcpos/2.0
    #ysc2 = np.random.random(num_imgs)*sourcpos-sourcpos/2.0
    #q = np.random.random(num_imgs)*0.5+0.5
    #vd = np.random.random(num_imgs)*100.0+200.0
    #pha = np.random.random(num_imgs)*360.0
    #zl = 0.2
    #zs = 1.0

    ysc1 = [-1.0]
    ysc2 = [0.0]
    zl = 0.2     #zl is the redshift of the lens galaxy.
    zs = 1.0
    vd = [520]    #Velocity Dispersion.
    q  = [0.5999999999999]
    pha = [0.0]


    #for i in xrange(rank,num_imgs,size):
    #for i in xrange(rank,num_imgs,size):
    for i in xrange(num_imgs):
        i = 0
        single_run_test(i,ysc1[i],ysc2[i],q[i],vd[i],pha[i],zl,zs)
