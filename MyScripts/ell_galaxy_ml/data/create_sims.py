import numpy as np
import pylab as pl
import astropy.io.fits as pyfits
import subprocess as sp


def gaussian_func_2d(x,y,mux,muy,sigmax,sigmay,theta):
	af0 = 1.0/(sigmax*sigmay*2.0*np.pi)
	af1 = np.cos(theta)**2.0/(2.0*sigmax**2.0)+np.sin(theta)**2.0/(2.0*sigmay**2.0)
	af2 =-np.sin(2.0*theta)/(4.0*sigmax**2.0)+np.sin(2.0*theta)/(4.0*sigmay**2.0)
	af3 = np.sin(theta)**2.0/(2.0*sigmax**2.0)+np.cos(theta)**2.0/(2.0*sigmay**2.0)

	res = af0*np.exp(-(af1*(x-mux)**2.0+2.0*af2*(x-mux)*(y-muy)+af3*(y-muy)**2.0))
	return res


def creat_cats(nimages):

    sigma = np.random.random(nimages)*1.5+0.5
    ell = np.random.random(nimages)*0.5+0.5
    pha = np.random.random(nimages)*360.0

    for i in xrange(nimages):
        print sigma[i], ell[i], pha[i]

    return 0


def create_imgs(par_cat):

    sigma, ell, pha = np.loadtxt(par_cat, usecols=(0,1,2), unpack=True)

    bsz = 9.0       #arcsec
    nnn = 50
    dsx = bsz/nnn   #arcsec

    x1 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
    x2 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
    x1, x2 = np.meshgrid(x1,x2)

    for i in xrange(len(sigma)):
        filename = "./fitsImages/"+str(i)+".fits"
        sigmax = sigma[i]*np.sqrt(ell[i])
        sigmay = sigma[i]/np.sqrt(ell[i])
        images = gaussian_func_2d(x1, x2, 0.0, 0.0, sigmax, sigmay, pha[i])
        pyfits.writeto(filename, images, clobber=True)


if __name__ == '__main__':
    nimgs = 1000
    #creat_cats(nimgs)
    create_imgs("./parameters.dat")
    #pl.show()
