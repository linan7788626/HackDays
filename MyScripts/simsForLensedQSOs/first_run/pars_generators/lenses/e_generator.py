import numpy as np
import pylab as pl
import scipy.interpolate as si


def pdf_e(e_var):
    e_mean = 0.3
    e_sigma = 0.16

    res0 = (e_var-e_mean)**2/e_sigma**2
    res = np.exp(-0.5*res0)
    return res


def e_generator(num_e):
    nbins = 5000
    randomeValues = np.random.random_sample(num_e)

    e_min = 0.0
    e_max = 0.9
    x = np.linspace(e_min, e_max, nbins)
    pdf = pdf_e(x)
    cdf = np.cumsum(pdf/np.sum(pdf))
    y = sorted(cdf)

    npdf = np.sum(pdf)*(x[1]-x[0])

    fdis =si.interp1d(y, x, kind='linear',bounds_error=False, fill_value=0.0)
    res = fdis(randomeValues)

    pl.figure()
    pl.hist(res, bins=40,range=(0,1), normed=1)
    pl.plot(x, pdf/npdf,'k-')

    return res


def pdf_e_tom(e_var,v_var):
    # e_var = 1.0-q_var
    s = 0.38+v_var*5.7e-4
    res = e_var/s**2.0*np.exp(-e_var*e_var/(2.0*s*s))
    return res


def e_generator_tom(num_e, vd):
    nbins = 50000
    randomeValues = np.random.random_sample(num_e)

    e_min = 0.0
    e_max = 0.8
    x = np.linspace(e_min, e_max, nbins)
    pdf = pdf_e_tom(x,vd)
    cdf = np.cumsum(pdf/np.sum(pdf))
    y = sorted(cdf)

    npdf = np.sum(pdf)*(x[1]-x[0])

    fdis =si.interp1d(y, x, kind='linear',bounds_error=False, fill_value=0.0)

    res = fdis(randomeValues)

    pl.figure()
    pl.hist(res, bins=40,range=(0,1), normed=1)
    pl.plot(x, pdf/npdf,'k-')

    return res

if __name__ == '__main__':
    e_generator(100000)
    # e_generator_tom(100000, 220)
    pl.show()
