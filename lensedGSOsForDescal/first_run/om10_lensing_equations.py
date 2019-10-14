import numpy as np
import pylab as pl

import om10
import triangle_root_finding as trf
import lambda_e as lef

from astropy.cosmology import FlatLambdaCDM
ncosmo = FlatLambdaCDM(H0=71, Om0=0.264, Ob0=0.044792699861138666)
vc = 2.998e5 #km/s
G = 4.3011790220362e-09 # Mpc/h (Msun/h)^-1 (km/s)^2
apr =  206269.43        #1/1^{''}
#--------------------------------------------------------------------

def Dc(z):
    res = ncosmo.comoving_distance(z).value*ncosmo.h
    return res


def Dc2(z1,z2):

    Dcz1 = (ncosmo.comoving_distance(z1).value*ncosmo.h)
    Dcz2 = (ncosmo.comoving_distance(z2).value*ncosmo.h)
    res = Dcz2-Dcz1+1e-8
    return res

def re_sv(sigmav, z1, z2):
    res = 4.0*np.pi*(sigmav/vc)**2.0*Dc2(z1, z2)/Dc(z2)*apr
    return res
#--------------------------------------------------------------------

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


def kappa_sie(x0, y0, theta, ql, re, le, x, y):
    tr = np.pi * (theta / 180.0)  # + np.pi / 2.0

    cs = np.cos(tr)
    sn = np.sin(tr)

    sx = x - x0
    sy = y - y0

    sx_r = sx * cs + sy * sn
    sy_r = -sx * sn + sy * cs

    res = re*le/2.0*1.0/np.sqrt(sx_r*sx_r*ql + sy_r*sy_r/ql)
    return res


def alphas_sie(x0, y0, theta, ql, re, le, ext_shears, ext_angle, ext_kappa, x, y):  # SIE lens model
    tr = np.pi * (theta / 180.0)   + np.pi / 2.0

    sx = x - x0
    sy = y - y0

    cs = np.cos(tr)
    sn = np.sin(tr)

    sx_r = sx * cs + sy * sn
    sy_r = -sx * sn + sy * cs

    eql = np.sqrt(ql / (1.0 - ql**2.0))
    psi = np.sqrt(sx_r**2.0 * ql + sy_r**2.0 / ql)
    dx_tmp = (re * eql * np.arctan( sx_r / psi / eql))
    dy_tmp = (re * eql * np.arctanh(sy_r / psi / eql))

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
    return dx*le + dx2 + dx3, dy*le + dy2 + dy3


def gauss_2d(x, y, xc, yc, sig):
    res0 = ((x-xc)**2+(y-yc)**2)/sig**2
    res = np.exp(-0.5*res0)
    return res


if __name__ == '__main__':
    print lef.lambda_e_tot(0.09182339)
    bsz = 4.0
    nnn = 512
    dsx = bsz/nnn

    x1, x2 = make_r_coor(nnn,dsx)
#---------------------------------------------
    db = om10.DB(catalog="/Users/uranus/GitHub/OM10/data/qso_mock.fits")
    # lid = 7176527
    # lid = 8519202
    # lid = 30184793
    # lid = 14864406
    # lid = 347938
    lid = 21393434
    lens = db.get_lens(lid)

    om10.plot_lens(lens)

    xl1 = 0.0
    xl2 = 0.0
    vd = lens.VELDISP[0]    # needed from OM10
    ql  = 1.0 - lens.ELLIP[0]
    le = lef.lambda_e_tot(lens.ELLIP[0])
    ph= lens.PHIE[0]
    zl = lens.ZLENS[0]
    zs = lens.ZSRC[0]

    ex_shs = lens.GAMMA[0]
    ex_sha = lens.PHIG[0]

    ys1 = lens.XSRC[0]
    ys2 = lens.YSRC[0]

    ximgs = lens.XIMG[0]
    yimgs = lens.YIMG[0]

    re0 = re_sv(vd, zl, zs)

    print vd, zl ,zs

    kap = kappa_sie(0.0, 0.0, ph, ql, re0, le, x1, x2)

    al1, al2 = alphas_sie(0.0, 0.0, ph, ql, re0, le, ex_shs, ex_sha, 0.0, x1, x2)

    print re0, le, ql, ph, ex_shs, ex_sha

    yi1 = x1 - al1
    yi2 = x2 - al2

    xroot1, xroot2, nroots = trf.mapping_triangles(ys1,ys2,x1,x2,yi1,yi2)

    if (nroots > len(ximgs[np.nonzero(ximgs)])):
        xroot = np.sqrt(xroot1*xroot1 + xroot2*xroot2)
        idx = xroot == xroot.min()
        xroot1[idx] = 0.0
        xroot2[idx] = 0.0

    simg = gauss_2d(x1,  x2,  ys1, ys2, 0.05)
    limg = gauss_2d(yi1, yi2, ys1, ys2, 0.05)

    levels_kappa = [0.8, 0.9, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.8]

    levels = [0.45, 0.6, 0.75, 0.9, 1.0]

    pl.figure(figsize=(8,8))
    pl.contourf(x1, x2, limg, levels)
    pl.contour(x1, x2, np.log(kap), levels_kappa, colors=('k',))
    pl.contour(x1, x2, simg, levels)
    pl.plot(ys1, ys2, 'rx')
    pl.plot(ximgs[np.nonzero(ximgs)], yimgs[np.nonzero(yimgs)], 'bo')
    pl.plot(xroot1[np.nonzero(xroot1)], xroot2[np.nonzero(xroot2)], 'rx')
#---------------------------------------------
    pl.show()
