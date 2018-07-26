#!/usr/bin/env python
import numpy as np
from astropy.cosmology import Planck13 as p13
import catalogize_lensed_images as sersic
apr = 206269.43

def return_re_le(sv_input,sv_array,re_array,le_array):
    dsv = 20.0
    idx = (np.abs(sv_array-sv_input)<dsv)
    Reff_r = np.random.choice(re_array[idx])
    log10_Ie_r = np.random.choice(le_array[idx])
    return Reff_r, log10_Ie_r #kpc, Lsun/pc^2

def lens_img_cat(lens_cat):
    xlc1 = lens_cat[0]
    xlc2 = lens_cat[1]
    ql = lens_cat[2]
    sigmav = lens_cat[3]
    phl = lens_cat[4]
    zl = lens_cat[5]

    zl_array, sv_array, re_array, le_array=np.loadtxt("./apj512628t1_ascii.txt",comments='#', usecols=(0,1,2,3),unpack=True)

    Reff_r,log10_Ieff_r = return_re_le(sigmav,sv_array, re_array, le_array)
    #Dl_l = p13.luminosity_distance(zl)
    Da_l = p13.angular_diameter_distance(zl).value
    Reff_arc = Reff_r/Da_l*apr*0.001 # arcsec
    Ieff_r = 10.0**(log10_Ieff_r) #/3.61 # http://www.astro.spbu.ru/staff/resh/Lectures/lec2.pdf
    ndex = 4
    mag_tot = sersic.sersic_mag_in_Rcut(Ieff_r,Reff_r,ndex)

    return [xlc1, xlc2, mag_tot, Reff_arc/np.sqrt(ql), Reff_arc*np.sqrt(ql), phl, ndex, zl]
if __name__ == '__main__':
    lens_cat = [0.0, 0.0, 0.5, 280, 47, 0.5]
    limg_cat = lens_img_cat(lens_cat)

    print limg_cat
