import numpy as np
import astropy.io.fits as pyfits

#---------------------------
# 'ACS WFC1' I UKIRT
PHOTFLAM = 7.0723600E-20    # erg cm-2 s-1 A-1
ZP_AB = 25.936
F0_AB = 2250.               # Jy
# EXP_TIME = 4056.            # s
EXP_TIME = 1.               # s
Jansky = 3e-13              # erg cm-2 s-1 A-1


def totimg_to_mag(inputs):

    # res = -2.5*np.log10(np.sum(inputs)*PHOTFLAM/Jansky/EXP_TIME)+ZP_AB

    res = -2.5*np.log10(inputs*PHOTFLAM/Jansky/EXP_TIME/F0_AB)

    return res


def mag_to_totimg(inputs):

    # res = -2.5*np.log10(np.sum(inputs)*PHOTFLAM/Jansky/EXP_TIME)+ZP_AB

    res = 10.0**(-inputs/2.5)/(PHOTFLAM/Jansky/EXP_TIME/F0_AB)

    return res


def zero_factors(totimg, totmag):

    totimg = 10.0**(-totmag/2.5)/(PHOTFLAM/Jansky/EXP_TIME/F0_AB)

    res = 10.0**(-totmag/2.5)/totimg

    return res

def img2_to_mag2(img1, mag1, img2):

    fz = 10.0**(-mag1/2.5)/img1
    res = -2.5*np.log10(img2*fz)

    return res

if __name__ == '__main__':
    img = pyfits.getdata('../../lens_modeling_challenge/gals_sources/439.0_149.482739_1.889989_processed.fits')
    mag = mag_cat[0]
    zl = zl_cat[0]

    tot = np.sum(img)

    mag = totimg_to_mag(tot)
    print tot, mag
    tot2 = mag_to_totimg(mag)
    print tot2, mag
