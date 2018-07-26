import numpy as np
import pylab as pl
import astropy.io.fits as pyfits

import om10
db = om10.DB(catalog="/Users/uranus/GitHub/OM10/data/qso_mock.fits")

# import subprocess as sp

apr = 206269.43 # arcseconds per rad


def cross_check_positions(lgals_uID):
    # prefix,uniqueId,raPhoSim,decPhoSim,phoSimMagNorm,sedFilepath,redshift,shear1,shear2,kappa,raOffset,decOffset,spatialmodel,internalExtinctionModel,galacticExtinctionModel,galacticAv,galacticRv,twinkles_system,twinkles_img_num,lens_galaxy_uID
    agn_ra, agn_dec, agn_twinkles_id, agn_lens_galaxy_uID = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(2,3,17,19), unpack=True, ndmin=0)

    agn_lens_galaxy_uID = agn_lens_galaxy_uID.astype('int')
    idx_agn = agn_lens_galaxy_uID == lgals_uID

    agn_ra = agn_ra.astype('double')
    agn_dec = agn_dec.astype('double')
    agn_twinkles_id = agn_twinkles_id.astype('int')

    lagn_ra = agn_ra[idx_agn]
    lagn_dec = agn_dec[idx_agn]

#----------------------------------------------------------------------------
# Parameters of the lens, by matching twinklesid
#
    twinklesid = agn_twinkles_id[idx_agn][0]
    hdulist = pyfits.open('./twinkles_DESC_SLAC/twinkles_lenses_v2.fits')
    # print hdulist[1].header   # needed from OM10

    idx = hdulist[1].data['twinklesId'] == twinklesid

    ximg = hdulist[1].data['XIMG'][idx][0]
    yimg = hdulist[1].data['YIMG'][idx][0]

    # TTYPE8  = 'GAMMA'   TTYPE9  = 'PHIG    '
#----------------------------------------------------------------------------
    #prefix,uniqueId,raPhoSim,decPhoSim,phoSimMagNorm,sedFilepath,redshift,shear1,shear2,kappa,raOffset,decOffset,spatialmodel,majorAxis,minorAxis,positionAngle,sindex,internalExtinctionModel,internalAv,internalRv,galacticExtinctionModel,galacticAv,galacticRv
    lens_gals_uid, lens_gals_ra, lens_gals_dec = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_lens_galaxies_230.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(1,2,3), unpack=True, ndmin=0)

    lens_gals_uid = lens_gals_uid.astype('int')
    lens_gals_ra = lens_gals_ra.astype('double')
    lens_gals_dec = lens_gals_dec.astype('double')

    idx_lgals = lens_gals_uid == lgals_uID
#--------------------------------------------------------------------
    lagn_ra = agn_ra[idx_agn]
    lagn_dec = agn_dec[idx_agn]
    lgals_ra = lens_gals_ra[idx_lgals]
    lgals_dec = lens_gals_dec[idx_lgals]

    xn_lagn = (lagn_ra-lgals_ra)*3600.0
    yn_lagn = (lagn_dec-lgals_dec)*3600.0
    print ximg
    print xn_lagn
    print yimg
    print yn_lagn

    pl.figure(figsize=(8,8))
    pl.plot(ximg[np.nonzero(ximg)], \
            yimg[np.nonzero(yimg)], 'bx')
    pl.plot(xn_lagn, yn_lagn, 'ro')

    return 0


if __name__ == '__main__':
    LGALS_UID = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_lens_galaxies_230.txt", dtype="int", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(1,), unpack=True, ndmin=0)

    for i in LGALS_UID[:10]:
        print i
        cross_check_positions(i)

    pl.show()
