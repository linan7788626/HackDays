#LB Cluster galaxy image!

import sys
import os
import math
import numpy
import logging
import time
import galsim
import getopt
import numpy as np
import pyfits


#@profile
def main(argv):

    #to call python draw_cluster_galaxies.py -i cat.cat  -o /data/bleeml/whatever -f g
    #read in our catalog, the location for out outputs and the filter
    cat_file_name = ''
    output_file_base = ''
    filt = ''

    opts, args = getopt.getopt(sys.argv[1:],"i:o:f:")

    print opts, args

    for opt, arg in opts:
      if opt == '-h':
          print 'draw_cluster_galaxies.py -i cat.cat  -o /data/bleeml/whatever -f g'
          sys.exit()
      elif opt in ("-i", "--ifile"):
          cat_file_name = arg
      elif opt in ("-o", "--ofile"):
          output_file_base = arg
      elif opt in ("-f", "--filter"):
          filt = arg

    output_file_name = output_file_base + '_'+filt + '.fits'

    print 'catalog file is :', cat_file_name
    print 'Output file is :', output_file_name
    print 'Filter is :', filt


    #Things we might want to change:

    pixel_scale = 0.2               # arcsec / pixel  (size units in input catalog are pixels)  0.09 is what we want.
    #zero_pt = 32.0
    #zero_pt = 22.0
    zero_pt = 30
    pixels_box = 256
    #image_size = 10.0/(pixel_scale/60.0)  #5 arcmin radius
    #image_size = (pixel_scale*pixels_box)/60/(pixel_scale/60.0)
    image_size = pixels_box
    center = image_size/2.0
    image = galsim.ImageF(image_size, image_size)

    if filt == 'g':
        dfilt  = 0         #0=g-band, 1=r-band, 2=i-band, 3=z-band
    if filt == 'r':
        dfilt  = 1         #0=g-band, 1=r-band, 2=i-band, 3=z-band
    if filt == 'i':
        dfilt  = 2         #0=g-band, 1=r-band, 2=i-band, 3=z-band
    if filt == 'z':
        dfilt  = 3         #0=g-band, 1=r-band, 2=i-band, 3=z-band


    gsparams = galsim.GSParams(maximum_fft_size=8192)     # max FFT size

    # In non-script code, use getLogger(__name__) at module scope instead.
    logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
    logger = logging.getLogger("drawing_cluster_gals")

    logger.info(' - pixel scale = %.2f,',pixel_scale)
    logger.info(' - zero_point = %.2f,',zero_pt)
    logger.info(' - catalog = %s,',cat_file_name)
    logger.info(' - output = %s,',output_file_name)
    logger.info(' - filter = %s,',filt)
    logger.info(' - dfilt = %.2f,',dfilt)

    # Read in the input catalog
    cat = galsim.Catalog(cat_file_name)

    for k in xrange(cat.nobjects-1):  #let's make this work with 1 galaxy!
        print k
    #for k in xrange(1):  #let's make this work with 1 galaxy!
        # Galaxy is a bulge + disk with parameters taken from the catalog:
        bulge_frac = cat.getFloat(k+1, 7 + dfilt)   #k+1 b/c text column header
        disk_HLR = cat.getFloat(k+1,12)*1.67835  #I think this is Mike's convention, Double Check!!
        bulge_n = cat.getFloat(k+1,6)
        if bulge_n > 6.2:   #catching galsim limits!
            bulge_n = 6.2

        if bulge_n < 0.3:
            bulge_n = 0.3

        bulge_HLR =cat.getFloat(k+1,11)
        apparent_mag = cat.getFloat(k+1,2+dfilt)
        gal_flux = math.pow(2.512,zero_pt-apparent_mag)
        gal_ellip = cat.getFloat(k+1,14)
        gal_PA = cat.getFloat(k+1,13)  #need to double check this convention against what Nan is giving me
        gal_x = cat.getFloat(k+1,0)/pixel_scale     #shifts are given in pixels
        gal_y = cat.getFloat(k+1,1)/pixel_scale

        if (gal_x > pixels_box/2) | (gal_x < -pixels_box/2) | (gal_y > pixels_box/2) | (gal_y < -pixels_box/2):
            continue

        #logger.info('Starting galaxy using:')
        #logger.info('- x,y,amag,n,BT,HLR,DHLR,PA,e (%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f)',center+gal_x,center+gal_y,apparent_mag,bulge_n,bulge_frac,bulge_HLR,disk_HLR,gal_PA,gal_ellip)
        #logger.info('- flux (%.1f)',gal_flux)

        bulge = galsim.Sersic(n=bulge_n, half_light_radius=bulge_HLR, gsparams=gsparams)

        # Objects may be multiplied by a scalar (which means scaling the flux) and also
        # added to each other.
        if bulge_frac < 1:
            #logger.info('Disk + Bulge')
            disk = galsim.Sersic(1,half_light_radius=disk_HLR, gsparams=gsparams)  #only make a disk if we need to; disk Sersic =1
            gal = bulge_frac * bulge + (1-bulge_frac) * disk
        else:
            gal = bulge


        # Set the overall flux of the combined object.
        gal = gal.withFlux(gal_flux)

        # Shear the object (ellipticity and PA)

        gal_shape = galsim.Shear(q=gal_ellip, beta=gal_PA*galsim.degrees)
        gal = gal.shear(gal_shape)

        # The center of the object is normally placed at the center of the postage stamp image.
        # You can change that with shift:  Are these in pixels or arcsec?

        #        if bulge_n < 4:
        #     stamp = gal.drawImage(scale=pixel_scale)
        #else:
        #    print 'photon shooting!'
        #stamp = gal.drawImage(scale=pixel_scale,method='phot')
        stamp = gal.drawImage(scale=pixel_scale,method='no_pixel')

        stamp.setCenter(gal_x+center, gal_y+center)
        bounds = stamp.bounds & image.bounds
        image[bounds] += stamp[bounds]

        #logger.info('Made galaxy profile')


        # Add Poisson noise to the image:
        #        image.addNoise(galsim.PoissonNoise(rng, sky_level * pixel_scale**2))

        #        logger.info('Drew image for object at row %d in the input catalog'%k)


 # Write the image to a file
    # Note: if the file already exists, this will overwrite it.

    #image.write(output_file_name)

    #nx,ny = np.shape(image.array)
    #print nx,ny
    #gg = image.array[nx/2-2048:nx/2+2048,ny/2-2048:ny/2+2048]
    gg = image.array
    nx,ny = np.shape(gg)
    print nx,ny
    pyfits.writeto(output_file_name,gg,clobber=True)


if __name__ == "__main__":
    main(sys.argv)
