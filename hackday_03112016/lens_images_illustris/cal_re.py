import pyfits
import numpy as np
import pylab as pl
import mycosmology as mm

data = pyfits.getdata("./Dropbox (T)/hackday_03112016/lens_images_illustris/test_i_new.fits")
data_flat = data.flatten()
data_tot = np.sum(data)
data_eff = 0.5*data_tot
idx_max = data == np.max(data)

bsz = 4.22
nnn = 512
dsx = bsz/nnn
x1 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
x2 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
x1, x2 = np.meshgrid(x1, x2)

x1max = x1[idx_max]
x2max = x2[idx_max]

rr = np.sqrt((x1-x1max)**2.0+(x2-x2max)**2.0)
rrf = rr.flatten()
idx = np.argsort(rrf)

rrf_sorted = rrf[idx]
data_flat_sorted_with_rrf = data_flat[idx]

data_cum = np.cumsum(data_flat_sorted_with_rrf)
data_dif = np.abs(data_cum-data_eff)

idx_eff = data_dif == np.min(data_dif)

print "r_eff = ", rrf_sorted[idx_eff]

print "r_eff_ph = ", rrf_sorted[idx_eff]*mm.Da(0.56)/mm.apr*1e3
