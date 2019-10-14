#!/usr/bin/env python
import numpy as np
from astropy.cosmology import Planck13 as p13
import sersic
apr = 206269.43


def return_re_le(sv_input,sv_array,re_array,le_array):
    dsv = 20.0
    idx = (np.abs(sv_array-sv_input)<dsv)
    Reff_r = np.random.choice(re_array[idx])
    log10_Ie_r = np.random.choice(le_array[idx])
    return Reff_r, log10_Ie_r #kpc, log10(Lsun/pc^2)


def vd_matching(sv_input,sv_array,re_array,le_array):
    abs_ddv = np.abs(sv_array-sv_input)
    idx = abs_ddv==min(abs_ddv)
    Reff_r = re_array[idx].mean()
    log10_Iebar_r = le_array[idx].mean()
    Iebar_r = 10.0**(log10_Iebar_r)
    return Reff_r, Iebar_r #kpc, Lsun/pc^2


def lens_img_cat(lens_cat):
    xlc1 = lens_cat[0]
    xlc2 = lens_cat[1]
    ql = lens_cat[2]
    sigmav = lens_cat[3]
    phl = lens_cat[4]
    zl = lens_cat[5]
    ndex = 4

    zl_array, sv_array, re_array, le_array=np.loadtxt("./apj512628t1_ascii.txt",comments='#', usecols=(0,1,2,3),unpack=True)
    Reff_r,Iebar_r = vd_matching(sigmav,sv_array, re_array, le_array)
    Ieff_r = sersic.sersic_Iebar_to_Ieff(Iebar_r,Reff_r,ndex)
    mag_tot = sersic.sersic_mag_in_Rcut(Ieff_r,Reff_r*1e3,ndex,zl)

    Da_l = p13.angular_diameter_distance(zl).value
    Reff_arc = Reff_r*0.001/Da_l*apr

    return [xlc1, xlc2, mag_tot, Reff_arc, ql, phl, ndex, zl]


if __name__ == '__main__':
    vds, zds, es, extshrs, extshas, lphis = np.loadtxt("/Users/uranus/Desktop/sl_sims_for_Joe/Lens-Population/pars_generators/lenses/lens_pop.cat", dtype='float', comments="#", skiprows=1, usecols=(0,1,2,3,4,5), unpack=True)

    print "### mag_tot, Reff, ql, lphi, ndex"

    for i in xrange(len(vds)):
        vd = vds[i]
        ql = 1.0-es[i]
        lp = lphis[i]
        zd = zds[i]

        lens_cat = [0.0, 0.0, ql, vd, lp, zd]
        limg_cat = lens_img_cat(lens_cat)

        print limg_cat[2], limg_cat[3], limg_cat[4], limg_cat[5], limg_cat[6]

    # header = "### mag_tot, Reff, ql, lphi, ndex"
    # np.savetxt("./lgal_pop.cat", np.transpose((), fmt='%.18e', delimiter=' ', newline='\n', header=header, footer='', comments='# ')

