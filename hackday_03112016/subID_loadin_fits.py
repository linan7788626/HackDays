#!/usr/bin/env python
"""  A very simple examle script to download, open, and plot an image
	 from the Illustris Image host site.

	 Author:  P. Torrey (ptorrey@mit.edu) 1/24/15

	 Update:  P. Torrey (ptorrey@mit.edu) 6/21/15
		w/ suggestions from Geferson Lucatelli and Connor Bottrell
"""
import numpy as np
import os
import sunpy.sunpy__load as sunpy__load		#
#import sunpy.sunpy__plot as sunpy__plot
import sunpy.sunpy__synthetic_image as ssi
import illustris_python as il
import astropy.io.fits as pyfits
import cosmocalc as mm
import pylab as pl
basePath="../Illustris-1"


def loadin_central_fits(haloID,zl,bsz,nnn,band='ACS_F606_NEW.res'):
	GroupFirstSub = il.groupcat.loadHalos(basePath,135,fields=['GroupFirstSub'])
	sub_galnrs = GroupFirstSub[haloID]

	my_api = "9dc591936258a12cb7192ebf79450272" # Input your own api for Illustris database
	dl_base = 'http://www.illustris-project.org'

	cmd = 'wget --content-disposition --header="API-Key: '+my_api+'" "'+dl_base+ \
	'/api/Illustris-1/snapshots/135/subhalos/'+str(sub_galnrs)+  \
	'/stellar_mocks/broadband.fits"'

	if( not (os.path.isfile("./broadband_"+str(sub_galnrs)+".fits")) ):
		os.system(cmd)
	filename = "./broadband_"+str(sub_galnrs)+".fits"

	hdulist = pyfits.open(filename)
	nx1 = hdulist[3].header['NAXIS1']
	nx2 = hdulist[3].header['NAXIS2']
	Dcmr_zl0 = mm.cosmocalc(0.01, H0=70.4, WM=0.2726, WV=0.7274)['DL_Mpc']
	Dcmr_zl1 = mm.cosmocalc(zl, H0=70.4, WM=0.2726, WV=0.7274)['DL_Mpc']
	bsz_init = np.rad2deg(hdulist[3].header['FOV'])*3600.0*Dcmr_zl0/Dcmr_zl1 # arcsec
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

	return_image_new = ssi.congrid(return_image_new,(nnn,nnn))
	return return_image_new
