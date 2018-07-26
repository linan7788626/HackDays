import numpy as np
import om10
import alens_arr as aa
import astropy.io.fits as pyfits
import cosmocalc as mm
import illustris_python as il
import subID_loadin_fits as slf
basePath = '../Illustris-1'
#---------------------------------------------------------------------------------
def get_subhalo_IDs(om10_id):
	db = om10.DB(catalog="/global/homes/l/linan/.local/lib/python2.7/site-packages/OM10/data/qso_mock.fits")

	lens = db.get_lens(om10_id)
	zl = lens.ZLENS[0]
	zs = lens.ZSRC[0]
	qe = 1.0/(1.0 - lens.ELLIP[0])
	phi= lens.PHIE[0]
	vd0 = lens.VELDISP[0]	 # needed from OM10
	mtotal0 = aa.m200_sv(vd0,zl)   # needed from OM10
	#----------------------------------------------------------------------
	# index, snapid,haloid,subhid,ag_stellar,ag_total,eq_stellar,eq_total,ms_stellar,ms_total,vd_stellar
	#
	snapid,haloid,subhid,ag_total,eq_stellar,ms_total,vd_stellar = np.loadtxt("./il1_results.dat",usecols=(1,2,3,5,6,9,10),unpack=True)

	nhalos = len(subhid)

	ksi_om10 = (ag_total/phi-1.0)**2.0/nhalos*0.0+(eq_stellar/qe-1.0)**2.0/nhalos*0.4 \
			 + (ms_total/mtotal0-1.0)**2.0/nhalos*0.01+(vd_stellar*(1.0+zl)/vd0-1.0)**2.0/nhalos*100.4

	idx = ksi_om10 == np.min(ksi_om10)

	snapid_n = int(snapid[idx])
	haloid_n = int(haloid[idx])
	subhid_n = int(subhid[idx])

	return snapid_n,haloid_n,subhid_n,zl

if __name__ == '__main__':

	om10_id = 7176527
	snapID, haloID, subhID, ZL = get_subhalo_IDs(om10_id)
	print snapID,haloID,subhID
	Dcmr_zl = mm.cosmocalc(ZL, H0=70.4, WM=0.2726, WV=0.7274)['DL_Mpc']
	bsz = 30.0/1000.0/Dcmr_zl*mm.arcsec_per_rad #arcsec
	nnn = 512

	lens_img = slf.loadin_central_fits(haloID,ZL,bsz,nnn,"g_SDSS.res")
	pyfits.writeto("test_g.fits",lens_img,clobber=True)
	lens_img = slf.loadin_central_fits(haloID,ZL,bsz,nnn,"r_SDSS.res")
	pyfits.writeto("test_r.fits",lens_img,clobber=True)
	lens_img = slf.loadin_central_fits(haloID,ZL,bsz,nnn,"i_SDSS.res")
	pyfits.writeto("test_i.fits",lens_img,clobber=True)
