import sys
del sys.path[2]
import om10
db = om10.DB(catalog="./qso_mock.fits")
id = 7176527
lens = db.get_lens(id)
print db.lenses.names

print "id =", lens.LENSID[0]
print "vd =", lens.VELDISP[0]    # needed from OM10
print "xi =", lens.XIMG[0]
print "yi =", lens.YIMG[0]
print "mui=", lens.MAG[0]
print "nim=", lens.NIMG[0]
print "md =", lens.APMAG_I[0]
print "ms =", lens.MAGI_IN[0]
print "xs =", lens.XSRC[0]
print "ys =", lens.YSRC[0]
print "xd =", 0.0
print "yd =", 0.0
print "zd =", lens.ZLENS[0]
print "zs =", lens.ZSRC[0]
print "q  =", 1.0 - lens.ELLIP[0]
print "phi=", lens.PHIE[0]
print "reff=", lens.REFF_T[0]

om10.plot_lens(lens)
#lenses = db.select_random(maglim=23.3,area=20000.0,IQ=0.7)
