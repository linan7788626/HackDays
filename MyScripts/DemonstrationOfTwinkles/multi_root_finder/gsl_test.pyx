cdef extern from "gsl/gsl_sf.h":
    double gsl_sf_gamma(double x)
    double GSL_SF_GAMMA_XMAX

def integrate_gamma(double a, double b,
                    int n=10000):
    if (min(a, b) <= 0 or
        max(a, b) >= GSL_SF_GAMMA_XMAX):
        raise ValueError(’Limits out ’
          ’of range (0, \%f)’ %
          GSL_SF_GAMMA_XMAX)
    cdef int i
    cdef double dx = (b - a) / n, result = 0
    for i in range(n):
        result += gsl_sf_gamma(a + i * dx) * dx
    return result
