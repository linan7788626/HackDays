import numpy as np
import scipy.special as ss
import scipy.integrate as sint

F = ss.hyp2f1

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
def Intergral_ud(llim,ulim,f,nbins=100000):

    soften = 1e-6
    #--------------------------------------
    # Scipy
    #
    # res = sint.romberg(f, llim+soften, ulim-soften)
    res = sint.quad(f, llim+soften, ulim-soften)[0]

    # #--------------------------------------
    # # Simplest
    # #
    # x = np.linspace(llim+soften, ulim-soften, nbins)
    # xl = x[:-1]
    # xu = x[1:]
    # dx = x[1] - x[0]

    # res = np.sum((f(xl)+f(xu))*dx/2.0)

    # #--------------------------------------
    # # Simpson's Rule
    # #
    # x = np.linspace(llim+soften, ulim-soften, nbins)
    # xl = x[:-1]
    # xu = x[1:]
    # xm = 0.5*(x[:-1]+x[1:])
    # dh = (x[1] - x[0])/2.0

    # res = np.sum((f(xl)+4.0*f(xm)+f(xu))*dh/3.0)

    return res


def lambda_e_obl(e):
    f2 = 1.0-e

    def WX_obl_qt(t, f2):
        q3 = np.sqrt(1.0-(1.0-f2*f2)/t**2)
        q3_p = np.sqrt((1.0-f2**2)/t**2)
        res = q3/q3_p*np.arcsin(q3_p)
        return res

    def WZ_obl_qt(t, f2):
        q3 = np.sqrt(1.0-(1.0-f2*f2)/t**2)
        q3_p = np.sqrt((1.0-f2**2)/t**2)
        res = 0.5*q3/(q3_p*np.arctan(q3_p/q3)) \
            * Intergral_ud(0.0, np.pi,
                       lambda theta: (1.0/np.sin(theta)
                                      *(np.arctan(q3_p/q3)**2.0
                                        -np.arctan(q3_p*np.cos(theta)/q3)**2.0)))
        return res

    res = np.pi/4.0*f2**1.5*F(1.0, 0.5, 2.0, f2**2) \
        / Intergral_ud(np.sqrt(1.0-f2*f2), 1.0,
                       lambda t: ((WX_obl_qt(t, f2)-WZ_obl_qt(t, f2))*t*t + WZ_obl_qt(t, f2))/np.sqrt(1.0-t*t))

    return res


def lambda_e_prol(e):
    f2 = 1.0-e

    def WX_prol_qt(t, f2):
        q3 = np.sqrt(1.0+(1.0/f2**2-1.0)/t**2)
        q3_p = np.sqrt((1.0/f2**2-1.0)/t**2)
        res = 0.5*q3/q3_p*np.log((q3+q3_p)/(q3-q3_p))
        return res

    def WZ_prol_qt(t, f2):
        q3 = np.sqrt(1.0+(1.0/f2**2-1.0)/t**2)
        q3_p = np.sqrt((1.0/f2**2-1.0)/t**2)
        res = 0.25*q3/(q3_p*np.log((q3+q3_p)/(q3-q3_p))) \
                * Intergral_ud(0.0,np.pi, lambda theta: (1.0/np.sin(theta) \
                                         * (np.log((q3+q3_p)/(q3-q3_p))**2 \
                                         - np.log((q3+q3_p*np.cos(theta)) \
                                                  /(q3-q3_p*np.cos(theta)))**2)))
        return res


    f2_p = np.sqrt(1.0-f2*f2)
    q3_max = 3.46717
    t_min = f2_p/f2/np.sqrt(q3_max**2.0-1.0)

    if t_min >= 1.0:
        t_min = 1.0/t_min

    res = Intergral_ud(t_min,1.0, \
                       lambda t: (np.sqrt((t*t*f2*f2+f2_p*f2_p)/(f2*(1-t*t)))/t)) \
                        / Intergral_ud(t_min,1.0, \
                       lambda t: (((WX_prol_qt(t, f2)-WZ_prol_qt(t, f2))*t*t \
                                   +WZ_prol_qt(t, f2))/np.sqrt(1-t*t)))
    return res


def lambda_e_tot(e, ef=0.5):
    if e < 1e-6:
        res = 1.0
    else:
        res = ef*lambda_e_prol(e) + (1.0 - ef)*lambda_e_obl(e)
    return res


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
