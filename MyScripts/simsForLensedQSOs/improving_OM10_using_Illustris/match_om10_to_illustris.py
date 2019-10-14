import sys
del sys.path[2]
import numpy as np
import pylab as pl
#import astropy.io.fits as pyfits
#import scipy.ndimage.filters as snf
#import scipy.signal as ss
import triangle_root_finding as trf
import alens_arr as aa
#import congrid
import scipy.interpolate as sci
#from scipy.ndimage.filters import gaussian_filter
import om10
import python_call_c_so as pcs
import illustris_python as il
import mycosmology as mm
basePath = '../Illustris-3'
#---------------------------------------------------------------------------------
def read_sdens_from_illustris(dm,xn1,xn2,xn3,bsz,nnn):
    dsx = bsz/nnn

    idx_ptcls_in1a = xn1>-bsz/2.0
    idx_ptcls_in1b = xn1<+bsz/2.0
    idx_ptcls_in2a = xn2>-bsz/2.0
    idx_ptcls_in2b = xn2<+bsz/2.0

    idx_ptcls_in = idx_ptcls_in1a&idx_ptcls_in1b&idx_ptcls_in2a&idx_ptcls_in2b

    xp1 = xn1[idx_ptcls_in]
    xp2 = xn2[idx_ptcls_in]
    xp3 = xn3[idx_ptcls_in]

    bzz = xp3.max()-xp3.min()
    dzz = bzz/nnn

    points = np.vstack([xp1,xp2,xp3]).T
    values = dm['SubfindDensity']*1e10
    values = values[idx_ptcls_in]

    grid_x, grid_y, grid_z = np.mgrid[(-bsz/2.0+dsx/2.0):(bsz/2.0-dsx/2.0):(nnn*1j),(-bsz/2.0+dsx/2.0):(bsz/2.0-dsx/2.0):(nnn*1j),xp3.min():xp3.max():(nnn*1j)]
    #xi1,xi2 = np.mgrid[0:1:(nnn*1j),0:1:(nnn*1j)]*bsz-bsz/2.0+dsx/2.0

    dens_3d = sci.griddata(points,values,(grid_x,grid_y,grid_z),method='linear',fill_value=0)
    sdens = np.sum(dens_3d,axis=2)*dzz*1e6
    return sdens
def make_r_coor(nc,dsx):

    bsz = nc*dsx
    x1 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0

    x2,x1 = np.meshgrid(x1,x2)
    return x1,x2
def make_c_coor(nc,dsx):

    bsz = nc*dsx
    x1,x2 = np.mgrid[0:(bsz-dsx):nc*1j,0:(bsz-dsx):nc*1j]-bsz/2.0+dsx/2.0
    return x1,x2

#--------------------------------------------------------------------
def point_ray_tracing(xi1,xi2,ai1,ai2,mu,spar):
    g_limage = xi1*0.0

    xroot1,xroot2,nroots = trf.roots_zeros(xi1,xi2,ai1,ai2,spar[0],spar[1])

    idr1 = ((np.array(xroot1)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    idr2 = ((np.array(xroot2)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    g_limage[idr1,idr2] = g_limage[idr1,idr2] + spar[2]*np.abs(mua[idr1,idr2])

    return g_limage

def calculate_new_ys(ximg1,ximg2,alpha1,alpha2,dsi):

    aimg1 = pcs.call_inverse_cic_single(alpha1,ximg1,ximg2,dsi)
    aimg2 = pcs.call_inverse_cic_single(alpha2,ximg1,ximg2,dsi)

    yimg1 = ximg1-aimg1
    yimg2 = ximg2-aimg2

    ysf1 = np.mean(yimg1)
    ysf2 = np.mean(yimg2)

    return ysf1,ysf2
#------------------------------------------------------------------------------
def get_halo_IDs(om10_id):
    db = om10.DB(catalog="./qso_mock.fits")

    lens = db.get_lens(om10_id)
    #nim= lens.NIMG[0]
    #md = lens.APMAG_I[0]
    #ms = lens.MAGI_IN[0]
    ys1 = lens.XSRC[0]
    ys2 = lens.YSRC[0]
    #xl = 0.0
    #yl = 0.0
    zl = lens.ZLENS[0]
    zs = lens.ZSRC[0]
    qe = 1.0/(1.0 - lens.ELLIP[0])
    phi= lens.PHIE[0]
    vd0 = lens.VELDISP[0]    # needed from OM10
    mtotal0 = aa.m200_sv(vd0,zl)   # needed from OM10
    print "Input: ", phi,qe,mtotal0,vd0
    #----------------------------------------------------------------------
    # index, snapid,haloid,subhid,ag_stellar,ag_total,eq_stellar,eq_total,ms_stellar,ms_total,vd_stellar
    #
    snapid,haloid,subhid,ag_total,eq_stellar,ms_total,vd_stellar = np.loadtxt("../produce_e_catalogs/results.dat",usecols=(1,2,3,5,6,9,10),unpack=True)

    nhalos = len(subhid)

    ksi_om10 = (ag_total/phi-1.0)**2.0/nhalos*0.0+(eq_stellar/qe-1.0)**2.0/nhalos*0.4 \
             + (ms_total/mtotal0-1.0)**2.0/nhalos*0.01+(vd_stellar*(1.0+zl)/vd0-1.0)**2.0/nhalos*100.4

    idx = ksi_om10 == np.min(ksi_om10)

    print "Matched: ",ag_total[idx],eq_stellar[idx],ms_total[idx],vd_stellar[idx]

    snapid_n = int(snapid[idx])
    haloid_n = int(haloid[idx])
    subhid_n = int(subhid[idx])

    lpar = [snapid_n,haloid_n,subhid_n,zl]
    spar = [ys1,ys2,lens.MAGI_IN[0],zs]
    #----------------------------------------------------------------------
    mui= lens.MAG[0]
    nim = lens.NIMG[0]
    ms = lens.MAGI_IN[0]
    muis = ms - 2.5*np.log10(np.abs(mui[:nim]))
    ximg1 = lens.XIMG[0][:nim]
    ximg2 = lens.YIMG[0][:nim]
    #print "om10.plot_lens: source magnitudes, magnification, image magnitudes:",ms,mui[:nim],mi
    ipar = [ximg1,ximg2,muis]
    return lpar,spar,ipar,lens

if __name__ == '__main__':

    om10_id = 7176527
    lpars,spars,ipars,lenses = get_halo_IDs(om10_id)
    #----------------------------------------------------------------------
    #om10.plot_lens(lenses)

    bsz = 50.0 # ckpc
    nnn = 512

    #GroupFirstSub = il.groupcat.loadHalos(basePath,snapid,fields=['GroupFirstSub'])
    stars = il.snapshot.loadSubhalo(basePath,lpars[0],lpars[2],'stars')
    nsp = stars['count']
    x1 = stars['Coordinates'][:,0]#/1e3/(1.0+zl) # Mpc/h, comoving
    x2 = stars['Coordinates'][:,1]#/1e3/(1.0+zl) # Mpc/h, comoving
    x3 = stars['Coordinates'][:,2]#/1e3/(1.0+zl) # Mpc/h, comoving


    ptnl_p = stars['Potential']
    idx_ptnl_p_min = ptnl_p == ptnl_p.min()

    xcp1 = x1[idx_ptnl_p_min]
    xcp2 = x2[idx_ptnl_p_min]
    xcp3 = x3[idx_ptnl_p_min]
    #sdens = read_sdens_from_illustris(stars,x1-xcp1,x2-xcp2,x3-xcp3,bsz,nnn)*(1.0+lpars[3])**2.0

    #pl.figure()
    #pl.plot(x1,x2,"r.")
    #pl.figure()
    #pl.plot(x1new,x2new,"b.")
    #pl.show()
    (x1new,x2new) = trf.xy_rotate(x1, x2, xcp1, xcp2, -14.0)
    sdens = read_sdens_from_illustris(stars,x1new,x2new,x3-xcp3,bsz,nnn)*(1.0+lpars[3])**2.0

    bsz = bsz/1e3/mm.Dc(lpars[3])*mm.apr #arcsec
    dsx = bsz/nnn #arcsec

    xi1,xi2  = make_r_coor(nnn,dsx)
    ai1,ai2 = pcs.call_cal_alphas0(sdens,nnn,bsz,lpars[3],spars[3])

    mua = pcs.call_alphas_to_mu(ai1,ai2,nnn,dsx)

    #om10.plot_lens(lenses)

    pl.figure()
    pl.contour(xi1,xi2,mua)
    pl.colorbar()

    pl.figure()
    pl.plot(spars[0],spars[1],'ko')
    pl.contour(xi1-ai1,xi2-ai2,mua)
    pl.colorbar()

    [ys1,ys2] = calculate_new_ys(ipars[0],ipars[1],ai1,ai2,dsx)
    xrt1,xrt2,nrts = trf.roots_zeros(xi1,xi2,ai1,ai2,ys1,ys2)
    #xrt1,xrt2,nrts = trf.roots_zeros(xi1,xi2,ai1,ai2,spars[0],spars[1])
    #mui = pcs.call_inverse_cic_single(mua,xrt1,xrt2,dsx)

    print nrts

    #print mui
    #print len(xrt1)
    #print xrt1

    mui= lenses.MAG[0]
    ms = lenses.MAGI_IN[0]
    lenses.NIMG[0] = nrts
    muis = ms - 2.5*np.log10(np.abs(mui[:nrts]))
    lenses.XIMG[0][:] = 0.0
    lenses.XIMG[0][:nrts] = xrt1
    lenses.YIMG[0][:] = 0.0
    lenses.YIMG[0][:nrts] = xrt2
    lenses.MAG[0][:nrts] = muis

    print nrts
    print muis


    om10.plot_lens(lenses)

    #pl.figure()
    #pl.plot(spars[0],spars[1],'ko')
    ##pl.plot(ipars[0],ipars[1],'gs')
    #pl.plot(xrt1,xrt2,'rs')
    #pl.show()
    #g_limage = point_ray_tracing(xi1,xi2,ai1,ai2,mua,spars)

    #output_filename = "./new_om10_imgs.fits"
    #pyfits.writeto(output_filename,g_limage,clobber=True)

    ##pl.show()
    ## Time Delay
    ##0.5*(1.0+0.17)/mm.vc*(mm.Da(0.17)*mm.Da(3.0)/mm.Da2(0.17,3.0))*mm.Mpc/1000.0*(20.0*2.0)/mm.apr**2.0/mm.yr * 365.0
