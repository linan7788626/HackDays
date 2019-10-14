import numpy as np
import pylab as pl
import scipy.special as ss
import scipy.interpolate as si

from astropy.cosmology import FlatLambdaCDM
ncosmo = FlatLambdaCDM(H0=71, Om0=0.264, Ob0=0.044792699861138666)
vc = 2.998e5 #km/s
G = 4.3011790220362e-09 # Mpc/h (Msun/h)^-1 (km/s)^2
apr =  206269.43        #1/1^{''}


def dn_dv(vd, zl=0.0):
    phi_star = 8e-3*(1.0+zl)**(0.229)
    vd_star = 161.0*(1.0+zl)**(-0.01)
    alpha = 2.23
    beta = 2.67

    res = phi_star*(vd/vd_star)**alpha*np.exp(-(vd/vd_star)**beta)*beta/(ss.gamma(alpha/beta)*vd)
    return res


def Nlens_dv(vd):
    vd_mean = 210.0
    vd_std = 54.0
    res0 = (vd-vd_mean)**2/(2.0*vd_std**2)
    res = np.exp(-res0)/np.sqrt(2.0*np.pi*vd_std**2.0)
    return res


def vd_generator(num_vd):
    nbins = 5000
    randomeValues = np.random.random_sample(num_vd)

    vd_min = 1.0
    vd_max = 400
    x = np.linspace(vd_min, vd_max, nbins)
    pdf = dn_dv(x)
    cdf = np.cumsum(pdf/np.sum(pdf))
    y = sorted(cdf)


    fdis =si.interp1d(y, x, kind='linear',bounds_error=False, fill_value=0.0)
    res = fdis(randomeValues)

    npdf = np.sum(pdf)*(x[1]-x[0])
    pl.figure()
    pl.hist(res, bins=40,range=(vd_min,vd_max), normed=1)
    pl.plot(x, pdf/npdf,'k-')

    return res


if __name__ == '__main__':

    # x = np.linspace(100.0,400.0,1000)
    # y = N_dv(x)

    # x1 = np.linspace(10.0,600.0,1000)
    # y1 = dn_dv(x1)

    # pl.figure(figsize=(8,5))
    # pl.plot(x, y, '-')
    # pl.figure(figsize=(8,5))
    # pl.plot(x1, y1, '-')
# #---------------------------------------------
    # x1 = np.linspace(10.0,600.0,1000)
    # y1 = dn_dv(x1, 0.0)
    # y2 = dn_dv(x1, 0.5)
    # y3 = dn_dv(x1, 1.5)

    # pl.figure(figsize=(8,5))
    # pl.plot(x1, y1, '-')
    # pl.plot(x1, y2, '-')
    # pl.plot(x1, y3, '-')
# #---------------------------------------------
    vd_generator(1000000)
# #---------------------------------------------
    # # vd = np.linspace(100.0,400.0,1000)
    # vd = 50.
    # zl = np.linspace(0.0,3.0,500)
    # dzl = zl[1] - zl[0]
    # dvol = ncosmo.differential_comoving_volume(zl)*dzl

    # N_dvol_dv = dn_dv(vd, zl)*dvol

    # # N_z = dn_dv(vd)*10.0*vol

    # pl.figure(figsize=(8,5))
    # pl.plot(zl, N_dvol_dv, '-')
# #---------------------------------------------
    pl.show()
