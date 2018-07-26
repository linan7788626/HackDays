import numpy as np
from astropy.cosmology import Planck13 as p13
import pyfits
apr = 206269.43

if __name__ == '__main__':
#----------------------------------------------------------------------------
    agn_id, agn_ra, agn_dec, agn_mag_norm, agn_redshift, agn_twinkles_id, agn_image_num, agn_lens_gal_id = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230.txt", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)
    agn_id = agn_id.astype('int')
    agn_twinkles_id = agn_twinkles_id.astype('int')
    agn_lens_gal_id = agn_lens_gal_id.astype('int')
#----------------------------------------------------------------------------
    hdulist = pyfits.open('./twinkles_DESC_SLAC/twinkles_lenses_v2.fits')
    twinklesid = agn_twinkles_id[:]

    xs = hdulist[1].data['XSRC'][twinklesid]   # needed from OM10
    ys = hdulist[1].data['YSRC'][twinklesid]
    zl = hdulist[1].data['ZLENS'][twinklesid]

#----------------------------------------------------------------------------
    host_id, host_ra, host_dec, host_mag_norm, host_redshift, host_major_axis, host_minor_axis = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_lens_galaxies_230.txt", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)
    host_id = host_id.astype('int')



    # Da_s = p13.angular_diameter_distance(agn_redshift).value
    # Da_l = p13.angular_diameter_distance(zl).value

    Dc_s = p13.comoving_distance(agn_redshift).value
    Dc_l = p13.comoving_distance(zl).value

    idx = agn_image_num == 0

    host_nid = np.array(agn_twinkles_id)[idx]
    lang_nid = np.array(agn_image_num)[idx]
    host_nra = np.array(xs)[idx]
    host_ndec= np.array(ys)[idx]
    host_nmag= np.array(np.random.choice(host_mag_norm, len(agn_id), replace=True)
                        - 5.0*np.log10(Dc_l/Dc_s))[idx]
    host_nred= np.array(agn_redshift)[idx]
    host_nma = np.array(np.random.choice(host_major_axis, len(agn_id)) * Dc_l/Dc_s)[idx]
    host_nmb = np.array(np.random.choice(host_minor_axis, len(agn_id)) * Dc_l/Dc_s)[idx]

    np.savetxt("./twinkles_DESC_SLAC/sprinkled_host_galaxies_230.txt",
               np.transpose((host_nid, lang_nid, host_nra, host_ndec, host_nmag, host_nred, host_nma, host_nmb)),
               fmt="%i, %i, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f", header="twinkleid,ra,dec,mag_norm,redshift,major_axis,minor_axis")
