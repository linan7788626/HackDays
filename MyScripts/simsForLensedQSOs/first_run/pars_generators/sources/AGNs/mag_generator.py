import numpy as np
import pylab as pl
import scipy.interpolate as si

from astropy.cosmology import FlatLambdaCDM
ncosmo = FlatLambdaCDM(H0=71, Om0=0.264, Ob0=0.044792699861138666)
vc = 2.998e5 #km/s
G = 4.3011790220362e-09 # Mpc/h (Msun/h)^-1 (km/s)^2
apr =  206269.43        #1/1^{''}
#--------------------------------------------------------------------
def pdf_e(M_var, z):

    xi = 2.98
    ki = 4.05
    z_star = 1.60

    fz = np.exp(xi*z)*(1.0+np.exp(ki*z_star)) \
        /(np.exp(ki*z/2)+np.exp(ki*z_star/2))**2.0

    M_star = -20.90 + 5.0*np.log10(ncosmo.h) - 2.5 * np.log10(fz)

    if (z <= 3) :
        alpha = -3.31
    else:
        alpha = -2.58
    beta = -1.45

    Phi_star = 5.34e-6
    res = Phi_star/(10.0**(0.4*(alpha+1)*(M_var-M_star)) \
                    +10.0**(0.4*(beta+1)*(M_var-M_star)))
    return res


def mag_generator(num_mag, z):
    nbins = 5000
    randomeValues = np.random.random_sample(num_mag)

    mag_min = -29.0
    mag_max = -24.0
    x = np.linspace(mag_min, mag_max, nbins)
    pdf = pdf_e(x, z)
    cdf = np.cumsum(pdf/np.sum(pdf))
    y = sorted(cdf)


    fdis =si.interp1d(y, x, kind='linear',bounds_error=False, fill_value=0.0)

    res = fdis(randomeValues)

    # npdf = np.sum(pdf)*(x[1]-x[0])
    # pl.figure()
    # pl.hist(res, bins=40,range=(mag_min, mag_max), normed=1)
    # pl.plot(x[::-1], np.log10(pdf/npdf),'k-')

    return res


if __name__ == '__main__':
    # bsz = 4.0
    # nnn = 512
    # dsx = bsz/nnn

    # x1, x2 = make_r_coor(nnn,dsx)
# #---------------------------------------------
    # db = om10.DB(catalog="/Users/uranus/GitHub/OM10/data/qso_mock.fits")
    # lid = 7176527
    # # lid = 8519202
    # # lid = 30184793
    # # lid = 14864406
    # lens = db.get_lens(lid)

    # print lens

    # om10.plot_lens(lens)

    # xl1 = 0.0
    # xl2 = 0.0
    # vd = lens.VELDISP[0]    # needed from OM10
    # le = 1.0
    # ql  = 1.0 - lens.ELLIP[0]
    # ph= lens.PHIE[0]
    # zl = lens.ZLENS[0]
    # zs = lens.ZSRC[0]

    # ex_shs = lens.GAMMA[0]
    # ex_sha = lens.PHIG[0]

    # ys1 = lens.XSRC[0]
    # ys2 = lens.YSRC[0]

    # ximgs = lens.XIMG[0]
    # yimgs = lens.YIMG[0]

    # # print vars(lens)
    # # print ximgs
    # # print yimgs

    # re0 = re_sv(vd, zl, zs)

    # kap = kappa_sie(0.0, 0.0, ph, ql, re0, le, x1, x2)

    # al1, al2 = alphas_sie(0.0, 0.0, ph, ql, re0, le, ex_shs, ex_sha, 0.0, x1, x2)

    # yi1 = x1 - al1
    # yi2 = x2 - al2

    # # xroot1, xroot2, nroots = trf.roots_zeros(x1, x2, al1, al2, ys1, ys2)
    # xroot1, xroot2, nroots = trf.mapping_triangles(ys1,ys2,x1,x2,yi1,yi2)

    # if (nroots > len(ximgs[np.nonzero(ximgs)])):
        # xroot = np.sqrt(xroot1*xroot1 + xroot2*xroot2)
        # idx = xroot == xroot.min()
        # xroot1[idx] = 0.0
        # xroot2[idx] = 0.0

    # print xroot1
    # print xroot2

    # simg = gauss_2d(x1,  x2,  ys1, ys2, 0.05)
    # limg = gauss_2d(yi1, yi2, ys1, ys2, 0.05)

    # # al11, al12 = np.gradient(al1, dsx)
    # # al21, al22 = np.gradient(al2, dsx)

    # # kap_tmp = (al11+al22)*0.5
    # # sh1_tmp = (al11-al22)*0.5
    # # sh2_tmp = (al12+al21)*0.5

    # # mu = 1.0/((1.0-kap_tmp)**2.0-sh1_tmp*sh1_tmp-sh2_tmp*sh2_tmp)

    # levels_kappa = [0.8, 0.9, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.8]
    # # pl.figure()
    # # pl.contour(x1, x2, mu, colors=('g',))
    # # pl.contour(x1-al1, x2-al2, mu, colors=('g',))
    # # pl.contour(x1, x2, np.log10(kap), levels, colors=('k',))
    # # pl.contour(x1, x2, np.log10(kap_tmp), levels, colors=('r',))
    # # pl.colorbar()

    # levels = [0.45, 0.6, 0.75, 0.9, 1.0]

    # pl.figure(figsize=(8,8))
    # pl.contourf(x1, x2, limg, levels)
    # pl.contour(x1, x2, np.log(kap), levels_kappa, colors=('k',))
    # pl.contour(x1, x2, simg, levels)
    # pl.plot(ys1, ys2, 'rx')
    # pl.plot(ximgs[np.nonzero(ximgs)], yimgs[np.nonzero(yimgs)], 'bo')
    # pl.plot(xroot1[np.nonzero(xroot1)], xroot2[np.nonzero(xroot2)], 'rx')
# #---------------------------------------------
    # mag_generator(10000, 1.25)
# #---------------------------------------------
    x = np.linspace(-29,-24, 1000)
    y1 = pdf_e(x, 1.25)
    y2 = pdf_e(x, 2.01)
    y3 = pdf_e(x, 3.75)

    pl.figure(figsize=(8,8))
    pl.xlim(-24, -29)
    pl.ylim(10**-10, 10**-4.8)
    pl.semilogy(x, y1*ncosmo.h**3.0/10, '-')
    pl.semilogy(x, y2*ncosmo.h**3.0, '-')
    pl.semilogy(x, y3*ncosmo.h**3.0*10, '-')
    pl.show()
