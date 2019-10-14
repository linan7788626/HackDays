import numpy as np
from astropy.cosmology import Planck13 as p13
from astropy import constants as const
import astropy.io.fits as pyfits
# import scipy.special as spf

import om10_lensing_equations as ole
import triangle_root_finding as trf
import lambda_e as lef
import libv4_cv as lv4
import pylab as pl
import generate_lens_images as gli
import sersic

apr = 206269.43 # arcseconds per rad

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


@profile
def main():
    bsz = 9.0
    nnn = 900
    dsx = bsz/nnn
    xi1, xi2 = make_r_coor(nnn,dsx)

    xlc1 = 0.0          # x position of the lens, arcseconds
    xlc2 = 0.0          # y position of the lens, arcseconds

    vd = 320.0
    ql = 0.7          # axis ratio b/a
    lp = 36.0           # orientation, degree
    ext_shears = 0.02
    ext_angle = 79.0
    zl = 0.5
    zs = 2.0          # redshift of the source

    rle = ole.re_sv(vd,zl,zs)
    #----------------------------------------------------------------------
    le = lef.lambda_e_tot(1.0-ql)

    ai1, ai2 = ole.alphas_sie(xlc1, xlc2, lp, ql, rle, le, ext_shears, ext_angle, 0.0, xi1, xi2)

    yi1 = xi1-ai1*1.0
    yi2 = xi2-ai2*1.0

    al12, al11 = np.gradient(ai1, dsx)
    al22, al21 = np.gradient(ai2, dsx)

    mua = 1.0/(1.0-(al11+al22)+al11*al22-al12*al21)
    #----------------------------------------------------------------------
    ysc1 = 0.0        # x position of the source, arcseconds
    ysc2 = 0.0        # y position of the source, arcseconds
    dsi = 0.03
    #----------------------------------------------------------------------
    mag_tot_sgal = 23.0     # total magnitude of the source
    g_source = pyfits.getdata(
        "../gals_sources/439.0_149.482739_1.889989_processed.fits")
    g_source = np.array(g_source.T, dtype="<d")
    g_source[g_source <= 0.0001] = 1e-6
    g_source = g_source.copy(order='C')

    g_limage = lv4.call_ray_tracing(g_source, yi1, yi2, ysc1, ysc2, dsi)
    mag_lensed_gal = img2_to_mag2(np.sum(g_source)*dsi*dsi,
                                  mag_tot_sgal,
                                  np.sum(g_limage)*dsx*dsx)
    #----------------------------------------------------------------------
    idac = g_source == g_source.max()

    yac1 = (np.indices(np.shape(g_source))[0][idac]-np.shape(g_source)[0]/2.0+0.5)*dsx
    yac2 = (np.indices(np.shape(g_source))[1][idac]-np.shape(g_source)[1]/2.0+0.5)*dsx

    # pl.figure()
    # pl.contourf(g_source.T, levels=[0.0,0.015,0.03,0.045,0.06,0.075,0.09,0.105])
    # pl.plot(np.indices(np.shape(g_source))[0][idac], np.indices(np.shape(g_source))[1][idac], 'rx')
    # pl.colorbar()

    mag_tot_sagn = 22.5
    a_source = g_source*0.0
    a_source[int((yac1+bsz/2.0-dsx/2.0)/dsx), int((yac2+bsz/2.0-dsx/2.0)/dsx)] = 1.0
    # a_source = a_source*mag2_to_img2(np.sum(g_source)*dsi*dsi, mag_tot_sgal, mag_tot_sagn)

    agns_rescale = mag2_to_img2(np.sum(g_source)*dsi*dsi, mag_tot_sgal, mag_tot_sagn) \
                /(np.sum(a_source)*dsi*dsi)

    a_source = a_source*agns_rescale

    a_limage = g_limage*0.0

    xroot1, xroot2, nroots = trf.mapping_triangles(yac1,yac2,xi1,xi2,yi1,yi2)

    # if (nroots > len(ximgs[np.nonzero(ximgs)])):
        # xroot = np.sqrt(xroot1*xroot1 + xroot2*xroot2)
        # idx = xroot == xroot.min()
        # xroot1[idx] = 0.0
        # xroot2[idx] = 0.0

    index2 = ((xroot1+bsz/2.0)/dsx).astype('int')
    index1 = ((xroot2+bsz/2.0)/dsx).astype('int')

    a_limage[index1, index2] = a_source.max()*dsi*dsi/(dsx*dsx)
    a_limage = a_limage*np.abs(mua)

    mag_lensed_agn = img2_to_mag2(np.sum(a_source)*dsi*dsi,
                                  mag_tot_sagn,
                                  np.sum(a_limage)*dsx*dsx)
    # print mag_lensed_agn
    # pl.figure()
    # pl.contourf(xi1, xi2, a_limage)
    # pl.plot(xroot1, xroot2, 'rx')
    # pl.colorbar()

    # print mag_lensed_gal
    # pl.figure()
    # pl.contourf(xi1, xi2, g_limage)
    # pl.plot(xroot1, xroot2, 'rx')
    # pl.colorbar()

    f_limage = g_limage + a_limage

    mag_lensed_tot = img2_to_mag2(np.sum(g_source)*dsi*dsi,
                                  mag_tot_sgal,
                                  np.sum(f_limage)*dsx*dsx)

    # print mag_lensed_tot
    # pl.figure()
    # pl.contourf(np.log10(f_limage))
    # pl.colorbar()

    #limg_cat = [xlc1, xlc2, mag_tot, Reff_arc, ql, phl, ndex, zl]

    lens_cat = [xlc1, xlc2, ql, vd, lp, zl]
    limg_cat = gli.lens_img_cat(lens_cat)

    limg_hst = sersic.sersic_2d_norm(xi1,xi2,limg_cat[0],limg_cat[1],limg_cat[3],ql,lp,4.0)

    limg_rescale = mag2_to_img2(np.sum(g_source)*dsi*dsi, mag_tot_sgal, limg_cat[2]) \
                /(np.sum(limg_hst)*dsx*dsx)

    g_lens = limg_hst*limg_rescale


    fimage = g_limage + a_limage + g_lens

    mag_tot = img2_to_mag2(np.sum(g_source)*dsi*dsi,
                                  mag_tot_sgal,
                                  np.sum(fimage)*dsx*dsx)

    fits_out = "./" + str(mag_tot_sgal) + "_" \
                    + str(mag_tot_sagn) + "_" \
                    + str(limg_cat[2]) + "_" \
                    + str(mag_tot) + "_full.fits"

    pyfits.writeto(fits_out, fimage, overwrite=True)
    # print mag_tot
    # pl.figure()
    # pl.contourf(np.log10(fimage))
    # pl.colorbar()

if __name__ == '__main__':
    main()
    pl.show()
