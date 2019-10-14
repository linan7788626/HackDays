#!/usr/bin/env python
"""  A very simple examle script to download, open, and plot an image
     from the Illustris Image host site.

     Author:  P. Torrey (ptorrey@mit.edu) 1/24/15

     Update:  P. Torrey (ptorrey@mit.edu) 6/21/15
		w/ suggestions from Geferson Lucatelli and Connor Bottrell
"""

#import sys
#del sys.path[2]
import numpy as np

import os
import sunpy.sunpy__load as sunpy__load		#
#import sunpy.sunpy__plot as sunpy__plot
import sunpy.sunpy__synthetic_image as ssi
import astropy.io.fits as pyfits
import mycosmology as mm
import pylab as pl


def loadin_fits(sub_galnrs,zl,bsz,nnn,band='ACS_F606_NEW.res'):

    my_api = "9dc591936258a12cb7192ebf79450272"
    dl_base='http://www.illustris-project.org'

    #try:
        #catalog = np.loadtxt('directory_catalog_135.txt',
                #dtype={'names'  : ('subdirs', 'galaxy_numbers', 'galaxy_masses'),
                               #'formats': ('S3', 'i8', 'f8')})
    #except:
        #os.system("wget "+dl_base+"/files/directory_catalog_135.txt")
        #catalog = np.loadtxt('directory_catalog_135.txt',
                #dtype={'names'  : ('subdirs', 'galaxy_numbers', 'galaxy_masses'),
                           #'formats': ('S3', 'i8','f8')})

    #1 rad2 = 3282.8 deg2 = 4.25 x 10^10 arcsec2

    #all_subdirs = catalog['subdirs']
    #all_galnrs  = catalog['galaxy_numbers']

    #int_all_subdirs = all_subdirs.astype('int')

    #sub_galnrs = 127228
    #sub_galnrs = 326247
    #sub_galnrs = 1030
    #sub_galnrs = 19

    cmd = 'wget --content-disposition --header="API-Key: '+my_api+'" "'+dl_base+ \
    '/api/Illustris-1/snapshots/135/subhalos/'+str(sub_galnrs)+  \
    '/stellar_mocks/broadband.fits"'

    if( not (os.path.isfile("./broadband_"+str(sub_galnrs)+".fits")) ):
        os.system(cmd)
    filename = "./broadband_"+str(sub_galnrs)+".fits"

    hdulist = pyfits.open(filename)
    nx1 = hdulist[3].header['NAXIS1']
    nx2 = hdulist[3].header['NAXIS2']
    bsz_init = np.rad2deg(hdulist[3].header['FOV'])*3600.0*mm.Dc(0.01)/mm.Dc(zl) # arcsec
    dsx_init = bsz_init/nx1
    hdulist.close()

    return_image = sunpy__load.load_broadband_image(filename,band)
    nnx1 = int(bsz/dsx_init)
    nnx2 = int(bsz/dsx_init)
    return_image_new = np.zeros((nnx1,nnx2))
    if bsz >= bsz_init:
        return_image_new[nnx1/2-nx1/2:nnx1/2+nx1/2,nnx2/2-nx2/2:nnx2/2+nx2/2] = return_image
    else:
        return_image_new = return_image[nx1/2-nnx1/2:nx1/2+nnx1/2,nx2/2-nnx2/2:nx2/2+nnx2/2]
    #print nx1/2-int(bsz/2.0/dsx_init),nx1/2+int(bsz/2.0/dsx_init),nx2/2-int(bsz/2.0/dsx_init),nx2/2+int(bsz/2.0/dsx_init)
    #return_image_new = return_image

    return_image_new = ssi.congrid(return_image_new,(nnn,nnn))
    #pl.figure()
    #pl.contourf(np.log10(return_image_new))
    #pl.colorbar()

    return return_image_new

if __name__ == '__main__':
    #bsz = 500.0/1000.0/mm.Dc(0.5)*mm.apr
    haloID = 0
    zl = 0.5
    bsz = 50.0
    nnn = 512
    lens_img = loadin_fits(haloID,zl,bsz,nnn)
    pyfits.writeto("fof.fits",lens_img,clobber=True)
    pl.show()
