import numpy as np
import pylab as pl
import scipy.interpolate as si


def pdf_extshr(e_var):
    e_mean = 0.05
    e_sigma = 0.2

    res = 1.0/np.sqrt(2.0*np.pi)/e_sigma/e_var \
        * np.exp(-(np.log(e_var)-e_mean)**2.0/(2.0*e_sigma**2.0))
    return res


def extshr_generator(num_e):
    nbins = 5000
    randomeValues = np.random.random_sample(num_e)

    e_min = 0.0001
    e_max = 2
    x = np.linspace(e_min, e_max, nbins)
    pdf = pdf_extshr(x)
    cdf = np.cumsum(pdf/np.sum(pdf))
    y = sorted(cdf)

    npdf = np.sum(pdf)*(x[1]-x[0])

    fdis =si.interp1d(y, x, kind='linear',bounds_error=False, fill_value=0.0)

    res = fdis(randomeValues)

    pl.figure()
    pl.hist(res, bins=40, normed=1)
    pl.plot(x, pdf/npdf,'k-')

    return res


if __name__ == '__main__':
    extshr_generator(10000)
    pl.show()
