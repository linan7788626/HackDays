http://www.stsci.edu/~dcoe/trilogy/Download.html
#------------------------------------------------------------------
http://www.stsci.edu/~dcoe/trilogy/Intro.html
#------------------------------------------------------------------
Current requirements are Python with the following libraries:

numpy

pyfits

PIL

scipy



To run:

alias trilogy python ~/trilogy/trilogy.py (or wherever the file resides)

trilogy trilogy.in



Convert any color image to an RGB Fits image viewable in ds9: im2rgbfits.py

usage: python im2rgbfits.py myimage.png <out_RGB.fits> <-header head.fits>

By default, the output image is myimage_RGB.fits, unless an output name is specified.

The output image will contain the header of another input file (here, head.fits), if given.



Convert a grayscale image to a regular FITS image: imgray2fits.py

#--------------------------------------------------------------
Trilogy uses log scaling constrained at three points (hence the name "tri-log"-y).

These three points are currently zero, "the noise", and saturated values.

The functional form used to accomplish this is y = log10( k * (x - xo) + 1 ) / r.



"The noise" is currently determined as 1-sigma above the sigma-clipped mean.

The user inputs how luminous "the noise" should be in the output image.

The default value of noiselum 0.15 works well, but you may wish to tweak this as high as 0.30 or so.



The user then inputs how what percentage of the data should be allowed to saturate.

satpercent = 0.001 (0.001% = 1/100,000 of the data) usually works well.

satpercent = 0 usually works well too.



Scaling is determined by sampling some region of the image: by default 1000x1000 at the center of the image. The user may need to adjust this using the parameters samplesize, sampledx, and sampledy, where the latter are offsets from the center.Â  For extragalactic work, we recommend including your brightest galaxies in this region while avoiding bright stars (allowing the latter to saturate).



To combine multiple filters in one channel (R,G,B), Trilogy currently simply adds the data.



Further documentation TBW. In the meantime, please e-mail me at dcoe@stsci.edu with questions.
