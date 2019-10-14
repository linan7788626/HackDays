# luminosity function of galaxies at different redshift.
# grab from SAM
# or CatSim
# Components: Sersic(Reff, Index), Mag-girb, ellipticity, orientation,

import numpy as np


def phi_m(M_var): # B-band
    phi_star = 1.6e-2
    M_star = -19.7+5.*np.log10(5.)
    alpha = -1.07

    res = np.log(10.0)/2.5*phi_star*10**(0.4*(alpha+1.)*(M_var-M_star)) \
        * np.exp(-10.0*(0.4*(M_var-M_star)))

    return res
