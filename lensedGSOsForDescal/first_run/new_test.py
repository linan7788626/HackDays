import numpy as np
from astropy.cosmology import Planck13 as p13
from astropy import constants as const
import astropy.io.fits as pyfits
from scipy.ndimage.filters import gaussian_filter
import pylab as pl
import ctypes as ct

import om10_lensing_equations as ole
import triangle_root_finding as trf0
import root_finding as trf
import lambda_e as lef
import congrid

apr = 206269.43 # arcseconds per rad

#---------------------------------------------------------------------------------
sps = ct.CDLL("../libsphsdens.so")

sps.cal_sph_sdens_weight.argtypes =[np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                    np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                    np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                    np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                    ct.c_float,ct.c_long,ct.c_float,ct.c_long,ct.c_long, \
                                    ct.c_float,ct.c_float,ct.c_float, \
                                    np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                    np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                    np.ctypeslib.ndpointer(dtype = ct.c_float)]

sps.cal_sph_sdens_weight.restype  = ct.c_int

def call_sph_sdens_weight(x1,x2,x3,mpp,Bsz,Nc,Np):

    x1 = np.array(x1,dtype=ct.c_float)
    x2 = np.array(x2,dtype=ct.c_float)
    x3 = np.array(x3,dtype=ct.c_float)
    mpp = np.array(mpp,dtype=ct.c_float)
    dsx = ct.c_float(Bsz/Nc)
    Ngb = ct.c_long(32)
    xc1 = ct.c_float(0.0)
    xc2 = ct.c_float(0.0)
    xc3 = ct.c_float(0.0)
    posx1 = np.zeros((Nc,Nc),dtype=ct.c_float)
    posx2 = np.zeros((Nc,Nc),dtype=ct.c_float)
    sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    sps.cal_sph_sdens_weight(x1,x2,x3,mpp,ct.c_float(Bsz),ct.c_long(Nc),dsx,Ngb,ct.c_long(Np),xc1,xc2,xc3,posx1,posx2,sdens);
    return sdens

#--------------------------------------------------------------------
# Construct regular grids
#
def make_r_coor(nc,dx):

    bs = nc*dx
    x1 = np.linspace(0,bs-dx,nc)-bs/2.0+dx/2.0
    x2 = np.linspace(0,bs-dx,nc)-bs/2.0+dx/2.0

    x1,x2 = np.meshgrid(x1,x2)
    return x1,x2
#--------------------------------------------------------------------
# Calculate Einstein Radius according to Velocity Dispersion
#
def re_sv(sv,z1,z2):
    Da_s = p13.angular_diameter_distance(z2).value
    Da_ls = p13.angular_diameter_distance_z1z2(z1,z2).value
    res = 4.0*np.pi*(sv**2.0/(const.c.value/1e3)**2.0)*Da_ls/Da_s*apr
    return res
#--------------------------------------------------------------------
# image total flux to magnitude and the reverse
#
def img2_to_mag2(img1_tot_flux, mag1, img2_tot_flux):
    fz = 10.0**(-mag1/2.5)/img1_tot_flux
    mag2 = -2.5*np.log10(img2_tot_flux*fz)
    return mag2


def mag2_to_img2(img1_tot_flux, mag1, mag2):
    fz = 10.0**(-mag1/2.5)/img1_tot_flux
    img2_tot_flux = 10.0**(-mag2/2.5)/fz
    return img2_tot_flux


def srcs_locations(bsz, nnn, xi1, xi2, yi1, yi2, mua):
    idx_crs = (mua<=-0.0)

    yc1_sample = yi1[idx_crs]
    yc2_sample = yi2[idx_crs]

    # ysc1 = np.random.choice(yc1_sample)        # x position of the source, arcseconds
    # ysc2 = np.random.choice(yc2_sample)        # y position of the source, arcseconds

    #----------------------------------------------------------------------
    # pl.figure()
    # pl.contour(xi1, xi2, mua, colors=('r',))
    # pl.plot(xi1[idx_crs], xi2[idx_crs], 'b.')

    #----------------------------------------------------------------------
    # Surface Density Parameters
    #

    #---------------------------
    # Calculate Surface Density
    #
    yc1_sample = yi1[idx_crs]
    yc2_sample = yi2[idx_crs]
    yc3 = np.linspace(-1.0, 1.0, len(yc1_sample))
    mp_in = np.ones((len(yc1_sample)))

    sdens = call_sph_sdens_weight(yc1_sample,yc2_sample,yc3,mp_in,bsz,nnn,len(yc1_sample))

    sth = 1000.

    ysc1 = np.random.choice(xi1[sdens > sth])
    ysc2 = np.random.choice(xi2[sdens > sth])

    # pl.figure(figsize=(8, 8))
    # pl.contour(yi1, yi2, mua, colors=('g',))
    # pl.plot(ysc1, ysc2,'k.')

    #---------------------------
    # Plot the contours
    #
    # pl.figure(figsize=(8, 8))
    # pl.contourf(np.log10(sdens))
    # pl.colorbar()

    map_tmp = np.zeros((nnn, nnn))
    map_tmp[sdens> sth ] = 1.0

    # pl.figure(figsize=(8, 8))
    # pl.contourf(xi1, xi2, map_tmp)
    # pl.contour(yi1, yi2, mua)
    # pl.plot(xi1[sdens > sth], xi2[sdens > sth], 'b.')

    return ysc1, ysc2


def pix2mag(pix_g_tot, pix_r_tot, pix_z_tot):
    fz_g = -3.31043359712e-09
    fz_r = -2.52341970111e-08
    fz_z = 1.52963399024e-08
    mag_g = -2.5*np.log10(pix_g_tot*fz_g)
    mag_r = -2.5*np.log10(pix_r_tot*fz_r)
    mag_z = -2.5*np.log10(pix_z_tot*fz_z)
    return mag_g, mag_r, mag_z


def mag2pix(mag_g, mag_r, mag_z):
    fz_g = -3.31043359712e-09
    fz_r = -2.52341970111e-08
    fz_z = 1.52963399024e-08
    img_tot_g = 10.0**(-mag_g/2.5)/fz_g
    img_tot_r = 10.0**(-mag_r/2.5)/fz_r
    img_tot_z = 10.0**(-mag_z/2.5)/fz_z
    return img_tot_g, img_tot_r, img_tot_z


def create_source_images(ysc1, ysc2, bsz, nnn, mag_in):

    dsi = bsz/nnn
    a_source = np.zeros((nnn, nnn))
    a_source[int((ysc1+bsz/2.0-dsi/2.0)/dsi), int((ysc2+bsz/2.0-dsi/2.0)/dsi)] = 1.0

    agns_rescale_g, agns_rescale_r, agns_rescale_z = mag2pix(mag_in)/np.sum(a_source)*dsi*dsi
    a_source_g = a_source*agns_rescale_g
    a_source_r = a_source*agns_rescale_r
    a_source_z = a_source*agns_rescale_z

    return a_source_g, a_source_r, a_source_z


def create_lensed_QSOs(a_source, ysc1, ysc2, xi1, xi2, yi1, yi2, mua, bsz, nnn):
    dsx = bsz/nnn

    a_limage = np.zeros((nnn, nnn))

    xroot1, xroot2, nroots = trf0.mapping_triangles(ysc1,ysc2,xi1,xi2,yi1,yi2)

    index2 = ((xroot1+bsz/2.0-dsx/2.0)/dsx).astype('int')
    index1 = ((xroot2+bsz/2.0-dsx/2.0)/dsx).astype('int')

    a_limage[index1, index2] = a_source.max()
    a_limage = a_limage*np.abs(mua)

    return a_limage


def mag2sigmav(mag_g, mag_r, zz):
    '''
    Calculate lens velocity dispersion using g,r magnitudes and redshift
    Use absolute magnitudes

    Get the SDSS magnitudes from megacam mags:
        g_Mega = g_SDSS - 0.153 (g_SDSS - r_SDSS)
        r_Mega = r_SDSS - 0.024 (g_SDSS - r_SDSS)
        g_Mega-r_Mega=(g_SDSS-r_SDSS)*(1-0.153+0.024)
        (g_Mega-r_Mega)/0.871=(g_SDSS-r_SDSS)
        r_SDSS = r_Mega + 0.024*(g_Mega-r_Mega)/0.871
    '''

    Dl_l = p13.luminosity_distance(zz).value
    mag_r_abs=mag_r-5.0*np.log10(Dl_l)-25.0

    ## Get the SDSS magnitude
    mag_r_sdss=mag_r_abs+0.024*(mag_g-mag_r)/0.871;

    ## Transfer to r' = ^{0.1}r Hubble type E, Frei and Gunn
    mag_r_sdss=mag_r_sdss-0.11;

    ## Assume that the luminosity function evolves such that Mr* decline by 1.5
    ## magnitudes from 0.0 to 1.0, this is similar to the evolution in the B band
    ## found by Faber et al. 2007 from DEEP2 and COMBO-17. This is really an adhoc
    ## prescription but nothing better known to me right now.

    mag_r_star=(-20.44)+(zz-0.1)*1.5;
    LbyLstar=10.0**(-0.4*(mag_r_sdss-mag_r_star));

    ## Parker et al. 2007, Table 1 - Bright sample
    res = 142.0*LbyLstar**(1./3.);

    return res


def final_lensed_QSOs_images(bsz, nnn, xi1, xi2, yi1, yi2, mua,\
                             mag_g_sagn, mag_r_sagn, mag_z_sagn, \
                             ysc1, ysc2, dsn, \
                             psf_g_width, psf_r_width, psf_z_width):

    dsx = bsz/nnn
    nnew = np.round(bsz/dsn)
    #----------------------------------------------------------------------
    xroots = trf.call_mapping_triangles([ysc1,ysc2],xi1,xi2,yi1,yi2)
    xrts = xroots[xroots.nonzero()]

    xroot1 = xrts[::2]
    xroot2 = xrts[1::2]

    index2 = ((xroot1+bsz/2.0-dsx/2.0)/dsx).astype('int')
    index1 = ((xroot2+bsz/2.0-dsx/2.0)/dsx).astype('int')
    #----------------------------------------------------------------------
    a_source = np.zeros((nnn, nnn))
    a_source[int((ysc1+bsz/2.0-dsx/2.0)/dsx), int((ysc2+bsz/2.0-dsx/2.0)/dsx)] = 1.0
    #----------------------------------------------------------------------
    agns_rescale_g = mag2pix(mag_g_sagn)/(np.sum(a_source)*dsx*dsx)
    g_source = a_source*agns_rescale_g
    g_limage = np.zeros((nnn, nnn))
    g_limage[index1, index2] = g_source.max()
    g_limage = g_limage*np.abs(mua)
    mag_g_lagn = pix2mag(np.sum(g_limage))
    g_limage_psf = gaussian_filter(g_limage, int(psf_g_width/dsx))
    g_limage_rebin = congrid.congrid(g_limage_psf, [nnew, nnew])
    rescale = mag2pix(mag_g_lagn)/np.sum(g_limage_rebin)*dsn*dsn
    g_limage_rebin_rescale = rescale*g_limage_rebin
    #----------------------------------------------------------------------
    agns_rescale_r = mag2pix(mag_r_sagn)/(np.sum(a_source)*dsx*dsx)
    r_source = a_source*agns_rescale_r
    r_limage = np.zeros((nnn, nnn))
    r_limage[index1, index2] = r_source.max()
    r_limage = r_limage*np.abs(mua)
    mag_r_lagn = pix2mag(np.sum(r_limage))
    r_limage_psf = gaussian_filter(r_limage, int(psf_r_width/dsx))
    r_limage_rebin = congrid.congrid(r_limage_psf, [nnew, nnew])
    rescale = mag2pix(mag_r_lagn)/np.sum(r_limage_rebin)*dsn*dsn
    r_limage_rebin_rescale = rescale*r_limage_rebin
    #----------------------------------------------------------------------
    agns_rescale_z = mag2pix(mag_z_sagn)/(np.sum(a_source)*dsx*dsx)
    z_source = a_source*agns_rescale_z
    z_limage = np.zeros((nnn, nnn))
    z_limage[index1, index2] = z_source.max()
    z_limage = z_limage*np.abs(mua)
    mag_z_lagn = pix2mag(np.sum(z_limage))
    z_limage_psf = gaussian_filter(z_limage, int(psf_z_width/dsx))
    z_limage_rebin = congrid.congrid(z_limage_psf, [nnew, nnew])
    rescale = mag2pix(mag_z_lagn)/np.sum(z_limage_rebin)*dsn*dsn
    z_limage_rebin_rescale = rescale*z_limage_rebin

    # # -------------------------------------------------------------
    # fits_out = "./" + str(mag_tot_sagn) + "_" \
                    # + str(mag_tot_lagn) + "_full.fits"

    # pyfits.writeto(fits_out, a_limage_rebin_rescale, overwrite=True)
    # # -------------------------------------------------------------

    # xi1_c = congrid.congrid(xi1, [nnew, nnew])
    # xi2_c = congrid.congrid(xi2, [nnew, nnew])

    # pl.figure()
    # pl.imshow(g_limage_rebin_rescale)
    # pl.colorbar()
    # pl.plot(xroot1, xroot2, 'ro')

    # pl.figure()
    # pl.contourf(r_limage_rebin_rescale)
    # pl.colorbar()

    # pl.figure()
    # pl.contourf(z_limage_rebin_rescale)
    # pl.colorbar()

    return g_limage_rebin_rescale, mag_g_lagn, \
           r_limage_rebin_rescale, mag_r_lagn, \
           z_limage_rebin_rescale, mag_z_lagn


def main(i):
    dsn = 0.262
    nno = 48
    index = i+1

    with pyfits.open('../data/decals_qso/decals_dr3_gal_dev.fits') as hdulist_cat:
        mag_g = hdulist_cat[1].data[i][3]
        mag_r = hdulist_cat[1].data[i][4]
        mag_z = hdulist_cat[1].data[i][5]

    with pyfits.open('../data/decals_qso/cutout/fits/gal/'+str(index)+'.fits') as hdulist_img:
        img_g = hdulist_img[0].data[0]
        # img_g[img_g<0.0] = 0.0
        img_r = hdulist_img[0].data[1]
        # img_r[img_r<0.0] = 0.0
        img_z = hdulist_img[0].data[2]
        # img_z[img_z<0.0] = 0.0

    img_tot_g, img_tot_r, img_tot_z = mag2pix(mag_g, mag_r, mag_z)
    mag_g_t, mag_r_t, mag_z_t = pix2mag(np.sum(img_g)*dsn*dsn, np.sum(img_r)*dsn*dsn, np.sum(img_z)*dsn*dsn)

    # print mag_g, mag_r, mag_z
    # print mag_g_t, mag_r_t, mag_z_t

    # print img_tot_g, img_tot_r, img_tot_z
    # print np.sum(img_g)*dsn*dsn, np.sum(img_r)*dsn*dsn, np.sum(img_z)*dsn*dsn

    fz_g = 10.0**(-mag_g/2.5)/(np.sum(img_g)*dsn*dsn)
    fz_r = 10.0**(-mag_r/2.5)/(np.sum(img_r)*dsn*dsn)
    fz_z = 10.0**(-mag_z/2.5)/(np.sum(img_z)*dsn*dsn)

    if ((fz_g<0.0) | (fz_r<0.0) | (fz_z<0.0)):
        print index

    # print fz_g
    # print fz_r
    # print fz_z

    # bsz = 21.0
    # nnn = 700
    # dsx = bsz/nnn
    # xi1, xi2 = make_r_coor(nnn,dsx)
    # #----------------------------------------------------------------------
    # # index, ra, dec, mag_g, mag_r, mag_z, ell, theta, psf_g, psf_r, psf_z, zl
    # with pyfits.open('../data/decals_qso/decals_dr3_gal_dev.fits') as hdulist_gal:
        # data_gal = hdulist_gal[1].data

    # # index, ra, dec, mag_g, mag_r, mag_z, psf_g, psf_r, psf_z, zs
    # with pyfits.open('../data/decals_qso/decals_dr3_qso_psf.fits') as hdulist_qso:
        # data_qso = hdulist_qso[1].data

    # mag_g = data_gal[i][3]
    # mag_r = data_gal[i][4]
    # # mag_z = data_gal[i][5]

    # xlc1 = 0.0          # x position of the lens, arcseconds
    # xlc2 = 0.0          # y position of the lens, arcseconds
    # zl = data_gal[i][11]

    # vd = mag2sigmav(mag_g, mag_r, zl)
    # # print vd
    # ql = 1.0-data_gal[i][6]         # axis ratio b/a
    # # print ql
    # lp = data_gal[i][7]          # orientation, degree
    # ext_shears = 0.02
    # ext_angle = 79.0

    # zs = data_qso[i][9]        # redshift of the source
    # rle = ole.re_sv(vd,zl,zs)
    # le = lef.e2le(1.0-ql)

    # ai1, ai2 = ole.alphas_sie(xlc1, xlc2, lp, ql, rle, le, ext_shears, ext_angle, 0.0, xi1, xi2)

    # yi1 = xi1-ai1
    # yi2 = xi2-ai2

    # al12, al11 = np.gradient(ai1, dsx)
    # al22, al21 = np.gradient(ai2, dsx)

    # mua = 1.0/(1.0-(al11+al22)+al11*al22-al12*al21)
    # #----------------------------------------------------------------------
    # ysc1, ysc2 = srcs_locations(bsz, nnn, xi1, xi2, yi1, yi2, mua)

    # mag_g_sagn = data_qso[i][3]
    # mag_r_sagn = data_qso[i][4]
    # mag_z_sagn = data_qso[i][5]
    # #----------------------------------------------------------------------
    # psf_g_width = data_qso[i][6]
    # psf_r_width = data_qso[i][7]
    # psf_z_width = data_qso[i][8]

    # #----------------------------------------------------------------------

    # lagn_g, lagn_g_mag, \
    # lagn_r, lagn_r_mag, \
    # lagn_z, lagn_z_mag = final_lensed_QSOs_images(bsz, nnn, xi1, xi2, yi1, yi2, mua, \
                                                  # mag_g_sagn, mag_r_sagn, mag_z_sagn, \
                                                  # ysc1, ysc2, dsn, \
                                                  # psf_g_width, psf_r_width, psf_z_width)

    return 0


if __name__ == '__main__':
    for i in xrange(5000):
        main(i)
    # pl.show()
