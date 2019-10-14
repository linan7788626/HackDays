import numpy as np

def Reff_Mv(M_var, z):
    # log10(Reff/kpc) = (Mv/-19.5)**-0.22*((1.0+z)/5.0)**-1.2 + N(0.3)
    res = 10.0**((M_var/-19.5)**-0.22*((1.0+z)/5.0)**-1.2 + np.random.normal(0.0, 1.0, 0.3))
    return res
