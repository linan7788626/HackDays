#!/usr/bin/env python
import numpy as np
import pyfits
apr = 206269.43

#----------------------------------------------------------------------------
def filter_out_img1():
    agn_id, agn_ra, agn_dec, agn_mag, agn_sed, agn_redshift, agn_catsim_id, agn_twinkles_id, agn_img_num = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230_catids.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)

    agn_id = agn_id.astype('int')
    agn_catsim_id = agn_catsim_id.astype('int')
    agn_twinkles_id = agn_twinkles_id.astype('int')
    agn_img_num = agn_img_num.astype('int')
#----------------------------------------------------------------------------
# Parameters of the lens, by matching twinklesid
#
    twinklesid = agn_twinkles_id
    hdulist = pyfits.open('./twinkles_DESC_SLAC/twinkles_lenses_v2.fits')

    print hdulist[1].header

    ys1 = hdulist[1].data['XSRC'][twinklesid]
    ys2 = hdulist[1].data['YSRC'][twinklesid]

    img_num = 0
    magnification_limg = np.abs(hdulist[1].data['MAG'][twinklesid])[:,0]

    idx = agn_img_num==img_num
    agnId = agn_id[idx]
    agnRa = agn_ra[idx]
    agnDec = agn_dec[idx]
    agnMag = agn_mag[idx]
    agnSed = agn_sed[idx]
    agnRedshift = agn_redshift[idx]
    agnCatsim_id = agn_catsim_id[idx]
    agnTwinkles_id = agn_twinkles_id[idx]
    agnImg_num = agn_img_num[idx]
    agnYs1 = ys1[idx]
    agnYs2 = ys2[idx]
    agnImgMag = magnification_limg[idx]

    # print "# AGN_ID, AGN_RA, AGN_DEC, AGN_ABMAG, AGN_SED, AGN_RedShift, AGN_CatSim_ID, AGN_Twinkles_ID, AGN_IMG_NUM, AGN_SRC_y1, AGN_SRC_y2, AGN_IMG_MAG"
    # for i in xrange(len(agnId)):
        # print agnId[i], agnRa[i], agnDec[i], agnMag[i], agnSed[i], agnRedshift[i], agnCatsim_id[i], agnTwinkles_id[i], agnImg_num[i], agnYs1[i], agnYs2[i], agnImgMag[i]

        # mag_src = host_mag[idx_host_components][0]-2.5*np.log10(1.0/magnification_limg)

        # srcs_cat = {'ys1' : ys1, 'ys2' : ys2, 'mag_src' : mag_src,
                    # 'Reff_src' : Reff_src, 'qs' : qs, 'phs' : phs,
                    # 'ns' : ns, 'zs' : zs, 'sed_src' : sed_src,
                    # 'agnid' : agnid, 'componentsid' : i}

    return 0
if __name__ == '__main__':
    filter_out_img1()
