import numpy as np
# import astropy.io.fits as pyfits
# import pylab as pl
import python_call_c_so as pcs
import triangle_root_finding as trf
import libv4_cv as lv4

def deflection_nie(xi1, xi2, xc1, xc2, re, rc, ql, theta, ext_shears, ext_angle, ext_kappa):  # SIE lens model
    tr = np.pi * (theta / 180.0)  # + np.pi / 2.0

    sx = xi1 - xc1
    sy = xi2 - xc2

    cs = np.cos(tr)
    sn = np.sin(tr)

    sx_r = sx * cs + sy * sn
    sy_r = -sx * sn + sy * cs

    psi = np.sqrt(ql**2.0 * (rc**2.0 + sx_r**2.0) + sy_r**2.0)
    dx_tmp = (re * np.sqrt(ql) / np.sqrt(1.0 - ql**2.0)) * \
        np.arctan(np.sqrt(1.0 - ql**2.0) * sx_r / (psi + rc))
    dy_tmp = (re * np.sqrt(ql) / np.sqrt(1.0 - ql**2.0)) * \
        np.arctanh(np.sqrt(1.0 - ql**2.0) * sy_r / (psi + rc * ql**2.0))
    dx = dx_tmp * cs - dy_tmp * sn
    dy = dx_tmp * sn + dy_tmp * cs

    # external shear
    tr2 = np.pi * (ext_angle / 180.0)
    cs2 = np.cos(2.0 * tr2)
    sn2 = np.sin(2.0 * tr2)
    dx2 = ext_shears * (cs2 * sx + sn2 * sy)
    dy2 = ext_shears * (sn2 * sx - cs2 * sy)

    # external kappa
    dx3 = ext_kappa * sx
    dy3 = ext_kappa * sy
    return dx + dx2 + dx3, dy + dy2 + dy3


def make_r_coor(nc,dsx):
    bsz = nc*dsx
    x1 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2,x1 = np.meshgrid(x1,x2)
    return x1,x2


def make_c_coor(nc,dsx):
    bsz = nc*dsx
    x1,x2 = np.mgrid[0:(bsz-dsx):nc*1j,0:(bsz-dsx):nc*1j]-bsz/2.0+dsx/2.0
    return x1,x2

#--------------------------------------------------------------------
def stacking_points_and_galaxies(xroot1,xroot2,lroots,lensed_imgs,bsz,nnn):
    res = lensed_imgs*0.0
    dsx = bsz/nnn

    idr1 = ((np.array(xroot1)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    idr2 = ((np.array(xroot2)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    res[idr1,idr2] = lroots

    res = res + lensed_imgs

    return res
#------------------------------------------------------------------------------
def run_main(lens_cat, lgal_cat, sagn_cat, sgal_cat):

    # lens_cat = {}
    # lgal_cat = {}
    # sagn_cat = {}
    # sgal_cat = {}

    bsz = 10.0 # arcsec
    nnn = 512
    dsx = bsz/nnn #arcsec
    dsi = 0.03

    zl = lens_cat['zl']
    zs = sagn_cat['zs']

    xc1= lens_cat['xc1']
    xc2= lens_cat['xc2']
    re = re_sv(lens_cat['SigmaV'], zl, zs)
    rc = 0.0
    ql = lens_cat['ql']
    phl = lens_cat['phl']
    esh = lens_cat['ext_shear']
    ash = lens_cat['ext_angle']
    ekp = lens_cat['ext_kappa']

    xi1, xi2 = make_r_coor(nnn, dsx)
    ai1, ai2 = deflection_nie(xi1, xi2, xc1, xc2, re, rc, ql, phl, esh, ash, ekp)
    mua = pcs.call_alphas_to_mu(ai1, ai2, nnn, dsx)

    yi1 = xi1 - ai1
    yi2 = xi2 - ai2
    #----------------------------------------------------------------------
    # Lensed Point Sources
    #
    ys1 = sagn_cat['ys1']
    ys2 = sagn_cat['ys2']
    xroot1, xroot2, nroots = trf.roots_zeros(xi1, xi2, ai1, ai2, ys1, ys2)
    muii = pcs.call_inverse_cic_single(mua, xroot1, xroot2, dsx)
    #------------------------------------------------------------------------------
    # Lensed Galaxies
    #
    ysc1 = ys1
    ysc2 = ys2

    sgal_img = loadin_sgal(sgal_cat)
    lensed_sgal = lv4.call_ray_tracing(sgal_img,yi1,yi2,ysc1,ysc2,dsi)
    #lensed_img = pcs.call_ray_tracing_single(ai1,ai2,nnn,bsz,zl)
    #------------------------------------------------------------------------------
    # Lens Galaxies
    #
    lgal_img = loadin_lgal(lgal_cat)

    #------------------------------------------------------------------------------
    # Stacking All images
    #
    res = stacking_points_and_galaxies(xroot1,xroot2,lroots[:,70],final_img,bsz,nnn)
    return res

if __name__ == '__main__':
    run_main()
