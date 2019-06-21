#!/usr/bin/env python
import numpy as np
import astropy.io.fits as pyfits
apr = 206269.43


def create_inputs_for_ray_tracing_agnid(agnid):
    agn_id, agn_ra, agn_dec, agn_mag, agn_sed, agn_redshift, agn_catsim_id, agn_twinkles_id, agn_img_num = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230_catids.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)

    agn_id = agn_id.astype('int')
    idx_agn = agn_id == agnid

    agn_catsim_id = agn_catsim_id.astype('int')
    agn_twinkles_id = agn_twinkles_id.astype('int')
    agn_img_num = agn_img_num.astype('int')
#----------------------------------------------------------------------------
    host_id, host_ra, host_dec, host_mag, host_sed, host_redshift, host_spatialmodel, host_major_axis, host_minor_axis, host_positionAngle, host_sindex, host_catsim_id = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_hosts_230.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)

    host_id = host_id.astype('int')
    host_catsim_id = host_catsim_id.astype('int')
    host_mag = host_mag.astype('double')
    host_redshift = host_redshift.astype('double')
    host_major_axis = host_major_axis.astype('double')
    host_minor_axis = host_minor_axis.astype('double')
    host_positionAngle = host_positionAngle.astype('double')
    host_sindex = host_sindex.astype('int')
    host_catsim_id = host_catsim_id.astype('int')
#----------------------------------------------------------------------------
# Parameters of the lens, by matching twinklesid
#
    twinklesid = agn_twinkles_id[idx_agn][0]
    hdulist = pyfits.open('./twinkles_DESC_SLAC/twinkles_lenses_v2.fits')
    # lid = hdulist[1].data['LENSID'][twinklesid]
    xl1 = 0.0
    xl2 = 0.0
    vd = hdulist[1].data['VELDISP'][twinklesid]   # needed from OM10
    zd = hdulist[1].data['ZLENS'][twinklesid]
    ql  = 1.0 - hdulist[1].data['ELLIP'][twinklesid]
    phi= hdulist[1].data['PHIE'][twinklesid]

    ys1 = hdulist[1].data['XSRC'][twinklesid]
    ys2 = hdulist[1].data['YSRC'][twinklesid]

    lens_cat = {'xl1' : xl1, 'xl2' : xl2, 'ql' : ql, 'vd' : vd,
                'phi' : phi, 'zd' : zd, 'twinklesid' : twinklesid}
#----------------------------------------------------------------------------
# Parameters of the sources
#
    idx2 = agn_img_num==0
    idx = idx_agn&idx2

    catsim_id = agn_catsim_id[idx]
    idx_host = host_catsim_id == catsim_id
    hostid = host_id[idx_host]
    srcs_cats = []
    for i in xrange(len(hostid)):
        srcs_cat = {}

        idx_host_components =  host_id == hostid[i]

        mag_src = host_mag[idx_host_components][0]
        Reff_src = np.sqrt(host_major_axis[idx_host_components][0]
                           *host_minor_axis[idx_host_components][0])
        qs = host_minor_axis[idx_host_components][0]/host_major_axis[idx_host_components][0]
        phs = host_positionAngle[idx_host_components][0]
        ns = host_sindex[idx_host_components][0]
        zs = host_redshift[idx_host_components][0]
        sed_src = host_sed[idx_host_components][0]

        srcs_cat = {'ys1' : ys1, 'ys2' : ys2, 'mag_src' : mag_src,
                    'Reff_src' : Reff_src, 'qs' : qs, 'phs' : phs,
                    'ns' : ns, 'zs' : zs, 'sed_src' : sed_src,
                    'agnid' : agnid, 'componentsid' : i}

        srcs_cats.append(srcs_cat)

        # lensed_images, lensed_mag, lensed_sed = ray_tracing(lens_cat, srcs_cat)

        # the output should be an image with TwinkleID + SED + Mag + zl +
    return lens_cat, srcs_cats
if __name__ == '__main__':
    agnid_tmp = 794901004316
    lens_cats_tmp, srcs_cats_tmp = create_inputs_for_ray_tracing_agnid(agnid_tmp)
  #  print lens_cats_tmp, srcs_cats_tmp


def create_inputs_for_ray_tracing():
#----------------------------------------------------------------------------
    lens_catt, srcs_catt = create_inputs_for_ray_tracing_agnid(agnid_tmp)

    return lens_catt, srcs_catt

if __name__ == '__main__':
    lens_catt_tmp, srcs_catt_tmp = create_inputs_for_ray_tracing()

    phile = open("lens.cat","w")
   # print >> phile, lens_cattt.values()
    #print >> phile, srcs_cattt.values()

 #   print lens_catt_tmp.values()[0]
    aaa = srcs_catt_tmp[0]
    bbb = srcs_catt_tmp[1]
 #   print bbb.values()[1]

    print aaa

    print >> phile, lens_catt_tmp.values()[5], lens_catt_tmp.values()[6], lens_catt_tmp.values()[1], lens_catt_tmp.values()[3], lens_catt_tmp.values()[0], lens_catt_tmp.values()[2], aaa.values()[7], aaa.values()[6], aaa.values()[0], aaa.values()[4], aaa.values()[5], aaa.values()[2], aaa.values()[8], aaa.values()[3]
    print >> phile, lens_catt_tmp.values()[5], lens_catt_tmp.values()[6], lens_catt_tmp.values()[1], lens_catt_tmp.values()[3], lens_catt_tmp.values()[0], lens_catt_tmp.values()[2], bbb.values()[7], bbb.values()[6], bbb.values()[0], bbb.values()[4], bbb.values()[5], bbb.values()[2], bbb.values()[8], bbb.values()[3]
#----------------------------------------------------------------------------
