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

    print len(values)

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
def point_ray_tracing(xi1,xi2,ai1,ai2,mua,spar):
    g_limage = xi1*0.0
    dsx = xi1[1,1]-xi1[0,0]
    bsz = np.shape(xi1)[0]*dsx

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
def get_halo_cc(om10_id):
    db = om10.DB(catalog="./qso_mock.fits")

    lens = db.get_lens(om10_id)
    zl = lens.ZLENS[0]
    zs = lens.ZSRC[0]
    qe = 1.0/(1.0 - lens.ELLIP[0])
    phi= lens.PHIE[0]
    vd0 = lens.VELDISP[0]    # needed from OM10
    mtotal0 = aa.m200_sv(vd0,zl)   # needed from OM10
    #----------------------------------------------------------------------
    # index, snapid,haloid,subhid,ag_stellar,ag_total,eq_stellar,eq_total,ms_stellar,ms_total,vd_stellar
    #
    snapid,haloid,subhid,ag_total,eq_stellar,ms_total,vd_stellar = np.loadtxt("../produce_e_catalogs/results.dat",usecols=(1,2,3,5,6,9,10),unpack=True)

    print "Input: ", qe, phi, vd0, mtotal0

    nhalos = len(subhid)

    ksi_om10 = (ag_total/phi-1.0)**2.0/nhalos*0.0+(eq_stellar/qe-1.0)**2.0/nhalos*0.4 \
             + (ms_total/mtotal0-1.0)**2.0/nhalos*0.01+(vd_stellar*(1.0+zl)/vd0-1.0)**2.0/nhalos*100.4

    idx = ksi_om10 == np.min(ksi_om10)
    print "Matched: ", eq_stellar[idx], ag_total[idx], vd_stellar[idx], ms_total[idx]

    snapid_n = int(snapid[idx])
    #haloid_n = int(haloid[idx])
    subhid_n = int(subhid[idx])

    bsz = 30.0 # ckpc
    nnn = 512

    stars = il.snapshot.loadSubhalo(basePath,snapid_n,subhid_n,'stars')
    x1 = stars['Coordinates'][:,0]#/1e3/(1.0+zl) # Mpc/h, comoving
    x2 = stars['Coordinates'][:,1]#/1e3/(1.0+zl) # Mpc/h, comoving
    x3 = stars['Coordinates'][:,2]#/1e3/(1.0+zl) # Mpc/h, comoving


    ptnl_p = stars['Potential']
    idx_ptnl_p_min = ptnl_p == ptnl_p.min()

    xcp1 = x1[idx_ptnl_p_min]
    xcp2 = x2[idx_ptnl_p_min]
    xcp3 = x3[idx_ptnl_p_min]

    (x1new,x2new) = trf.xy_rotate(x1, x2, xcp1, xcp2, -14.0)
    sdens = read_sdens_from_illustris(stars,x1new,x2new,x3-xcp3,bsz/512*10.0,nnn)*(1.0+zl)**2.0

    kappa = sdens/mm.sigma_crit(zl,zs)

    bsz = bsz/1e3/mm.Dc(zl)*mm.apr #arcsec
    dsx = bsz/nnn #arcsec

    xi1,xi2  = make_r_coor(nnn,dsx)
    ai1,ai2 = pcs.call_cal_alphas0(sdens,nnn,bsz,zl,zs)

    mua = pcs.call_alphas_to_mu(ai1,ai2,nnn,dsx)
    return xi1,xi2,ai1,ai2,mua,kappa

def lensing_model_SIE(x1,x2,lpar):
    xc1 = lpar[0]      # x coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]      # y coordinate of the center of lens (in units of Einstein radius).
    q   = lpar[2]      # Axis ratio of lens.
    rc  = lpar[3]      # Core size of lens (in units of Einstein radius).
    re  = lpar[4]      # Einstein radius of lens (in units of arcsec).
    pha = lpar[5]      # Orientation of lens (in units of degree).

    phirad = np.deg2rad(pha)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1-xc1)*cosa+(x2-xc2)*sina
    xt2 = (x2-xc2)*cosa-(x1-xc1)*sina

    phi = np.sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc)
    sq = np.sqrt(1.0-q*q)
    pd1 = phi+rc/q
    pd2 = phi+rc*q
    fx1 = sq*xt1/pd1
    fx2 = sq*xt2/pd2
    qs = np.sqrt(q)

    a1 = qs/sq*np.arctan(fx1)
    a2 = qs/sq*np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 =-sina

    fx11 = xt11/pd1-xt1*(xt1*q*q*xt11+xt2*xt21)/(phi*pd1*pd1)
    fx22 = xt22/pd2-xt2*(xt1*q*q*xt12+xt2*xt22)/(phi*pd2*pd2)
    fx12 = xt12/pd1-xt1*(xt1*q*q*xt12+xt2*xt22)/(phi*pd1*pd1)
    fx21 = xt21/pd2-xt2*(xt1*q*q*xt11+xt2*xt21)/(phi*pd2*pd2)

    a11 = qs/(1.0+fx1*fx1)*fx11
    a22 = qs/(1.0-fx2*fx2)*fx22
    a12 = qs/(1.0+fx1*fx1)*fx12
    a21 = qs/(1.0-fx2*fx2)*fx21

    rea11 = (a11*cosa-a21*sina)*re
    rea22 = (a22*cosa+a12*sina)*re
    rea12 = (a12*cosa-a22*sina)*re
    rea21 = (a21*cosa+a11*sina)*re

    kappa = 0.5*(rea11+rea22)
    shear1 = 0.5*(rea12+rea21)
    shear2 = 0.5*(rea11-rea22)

    y11 = 1.0-rea11
    y22 = 1.0-rea22
    y12 = 0.0-rea12
    y21 = 0.0-rea21

    jacobian = y11*y22-y12*y21
    mu = 1.0/jacobian

    alpha1 = (a1*cosa-a2*sina)*re
    alpha2 = (a2*cosa+a1*sina)*re

    return alpha1,alpha2,kappa,shear1,shear2,mu

def get_om10_cc(xi1,xi2,om10_id):
    db = om10.DB(catalog="./qso_mock.fits")

    lens = db.get_lens(om10_id)
    zl = lens.ZLENS[0]
    zs = lens.ZSRC[0]
    qe = 1.0/(1.0 - lens.ELLIP[0])
    phi= lens.PHIE[0]
    vd0 = lens.VELDISP[0]    # needed from OM10

    re = aa.re_sv(vd0,zl,zs)

    lpars = np.zeros((6))

    lpars[0] = 0.0 # x coordinate of the center of lens (in units of Einstein radius).
    lpars[1] = 0.0 # y coordinate of the center of lens (in units of Einstein radius).
    lpars[2] = 1.0/qe     # Axis ratio of lens.
    lpars[3] = 0.0     # Core size of lens (in units of Einstein radius).
    lpars[4] = re     # Einstein radius of lens (in units of arcsec).
    lpars[5] = phi     # Orientation of lens (in units of degree).

    alpha1, alpha2, kappa, shears1,shear2,mu = lensing_model_SIE(xi1,xi2,lpars)
    return alpha1,alpha2,mu,kappa

if __name__ == '__main__':

    om10_id = 7176527
    xi1,xi2,ai1,ai2,mui,kappai = get_halo_cc(om10_id)
    ao1,ao2,muo,kappao = get_om10_cc(xi1,xi2,om10_id)
    kappah = np.fromfile("./il1_kappa.bin",dtype=np.double)
    kappah = kappah.reshape((512,512))

    dsx = xi1[1,1]-xi1[0,0]
    bsz = dsx*512
    xbin = dsx*np.linspace(0,256-1,256)+dsx*0.5

    pl.figure()
    pl.ylim(0.1,100.0)
    pl.xlabel("arcsec")
    pl.ylabel(r"$\kappa$")
    pl.loglog(xbin,kappai[255,256:],'r',label="Illustris-3")
    pl.loglog(xbin,kappah[255,256:],'b',label="Illustris-1")
    pl.loglog(xbin,kappao[255,256:],'k',label="OM10")
    pl.legend()

    #pl.figure(figsize=(8,8))
    #pl.xlim(-1,1)
    #pl.ylim(-1,1)
    #pl.xlabel("arcsec")
    #pl.ylabel("arcsec")
    #pl.contour(xi1,xi2,mui,colors=['r',])
    #pl.contour(xi1,xi2,muo,colors=['k',])

    #pl.figure(figsize=(8,8))
    #pl.xlim(-0.5,0.5)
    #pl.ylim(-0.5,0.5)
    #pl.xlabel("arcsec")
    #pl.ylabel("arcsec")
    #pl.contour(xi1-ai1,xi2-ai2,mui,colors=['r',])
    #pl.contour(xi1-ao1,xi2-ao2,muo,colors=['k',])
    pl.show()
