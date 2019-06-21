#!/usr/bin/env python
import numpy as np
from astropy.cosmology import Planck13 as p13
import catalogize_lensed_images as sersic
import pyfits
apr = 206269.43

def return_re_le(sv_input,sv_array,re_array,le_array):
    dsv = 20.0
    idx = (np.abs(sv_array-sv_input)<dsv)
    Reff_r = np.random.choice(re_array[idx])
    log10_Ie_r = np.random.choice(le_array[idx])
    return Reff_r, log10_Ie_r #kpc, log10(Lsun/pc^2)

def vd_matching(sv_input,sv_array,re_array,le_array):
    abs_ddv = np.abs(sv_array-sv_input)
    idx = abs_ddv==min(abs_ddv)
    Reff_r = re_array[idx].mean()
    log10_Iebar_r = le_array[idx].mean()
    Iebar_r = 10.0**(log10_Iebar_r)
    return Reff_r, Iebar_r #kpc, Lsun/pc^2

def lens_img_cat(lens_cat):
    xlc1 = lens_cat[0]
    xlc2 = lens_cat[1]
    ql = lens_cat[2]
    sigmav = lens_cat[3]
    phl = lens_cat[4]
    zl = lens_cat[5]
    ndex = 4

    zl_array, sv_array, re_array, le_array=np.loadtxt("./apj512628t1_ascii.txt",comments='#', usecols=(0,1,2,3),unpack=True)
    Reff_r,Iebar_r = vd_matching(sigmav,sv_array, re_array, le_array)
    Ieff_r = sersic.sersic_Iebar_to_Ieff(Iebar_r,Reff_r,ndex)
    mag_tot = sersic.sersic_mag_in_Rcut(Ieff_r,Reff_r*1e3,ndex,zl)

    Da_l = p13.angular_diameter_distance(zl).value
    Reff_arc = Reff_r*0.001/Da_l*apr
    return [xlc1, xlc2, mag_tot, Reff_arc/np.sqrt(ql), Reff_arc*np.sqrt(ql), phl, ndex, zl, lens_cat[6]]

if __name__ == '__main__':
    ind = 0
#----------------------------------------------------------------------------
    agn_id, agn_ra, agn_dec, agn_mag_norm, agn_redshift, agn_twinkles_id, agn_image_num, agn_lens_gal_id = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230.txt", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)
    agn_id = agn_id.astype('int')
    agn_twinkles_id = agn_twinkles_id.astype('int')
    agn_lens_gal_id = agn_lens_gal_id.astype('int')
# #----------------------------------------------------------------------------
    # host_id, host_ra, host_dec, host_mag_norm, host_redshift, host_major_axis, host_minor_axis = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_lens_galaxies_230.txt", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)
    # host_id = host_id.astype('int')
#----------------------------------------------------------------------------
    hdulist = pyfits.open('./twinkles_DESC_SLAC/twinkles_lenses_v2.fits')

    twinklesid = agn_twinkles_id[0]

    lid = hdulist[1].data['LENSID'][twinklesid]
    xs = hdulist[1].data['XSRC'][twinklesid]
    ys = hdulist[1].data['YSRC'][twinklesid]
    vd = hdulist[1].data['VELDISP'][twinklesid]   # needed from OM10
    zd = hdulist[1].data['ZLENS'][twinklesid]
    ql  = 1.0 - hdulist[1].data['ELLIP'][twinklesid]
    phi= hdulist[1].data['PHIE'][twinklesid]

    xl1 = 0.0
    xl2 = 0.0
    lens_cat = [xl1, xl2, vd, ql, phi, zd, twinklesid]

    print "Lens: [xl1, xl2, vd, ql, phl, zd, twinkle_id]"
    print lens_cat

    host_id, host_ra, host_dec,host_mag_norm, host_redshift, host_major_axis, host_minor_axis, host_position_angle, host_sindex, host_catsim_id = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_hosts_230.txt", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(0, 1, 2, 3, 5, 7,8,9,10,11), unpack=True, ndmin=0)

    mag_src = host_mag_norm[0]
    major = host_major_axis[0]
    minor = host_minor_axis[0]
    phs = host_position_angle[0]
    zs = host_redshift[0]
    ndex = host_sindex[0]

    host_cat = [xs, ys, mag_src, np.sqrt(major*minor), minor/major, phs, ndex, zs]
    print "Host: [ysc1, ysc2, mag_tot, Reff, qs, phs, ndex, zs]"
    print host_cat

#----------------------------------------------------------------------------
    agn_id, agn_ra, agn_dec, agn_mag_norm, agn_redshift, agn_twinkles_id, agn_image_num, agn_lens_gal_id = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230.txt", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)
    agn_id = agn_id.astype('int')
    agn_twinkles_id = agn_twinkles_id.astype('int')
    agn_lens_gal_id = agn_lens_gal_id.astype('int')
#----------------------------------------------------------------------------
    lgal_id, lgal_ra, lgal_dec, lgal_mag_norm, lgal_redshift, lgal_major_axis, lgal_minor_axis = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_lens_galaxies_230.txt", dtype="string", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)
    lgal_id = lgal_id.astype('int')


    xc1 = lgal_ra[0]
    xc2 = lgal_dec[0]
    mag_tot = lgal_mag_norm[0]
    Reff_a = lgal_major_axis[0]
    Reff_b = lgal_minor_axis[0]
    ndex = 4.0
    phl = phi
    zl = zd

    lgal_cat = [xc1, xc2, mag_tot, Reff_a, Reff_b, phl, ndex, zl]

#----------------------------------------------------------------------------
    print "Lgal: [xc1, xc2, mag_tot, Reff_a, Reff_b, Theta, ndex, zl]"
    print lgal_cat

#----------------------------------------------------------------------------
