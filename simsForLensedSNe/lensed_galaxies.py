#/usr/bin/env python
import sys
del sys.path[2]
import numpy as np
import astropy.io.fits as pyfits
import pylab as pl
import python_call_c_so as pcs
import illustris_python as il
import mycosmology as mm
import scipy.interpolate as sci
import subID_loadin_fits as slf
#import libv4_cv as lv4


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
def run_main():
    HaloID = 0
    SnapID = 135

    bsz = 500.0 # ckpc
    zl = 0.5
    zs0 = 10.0
    nnn= 512

    dm = il.snapshot.loadHalo(basePath,SnapID,HaloID,'dm')
    #npp = dm['count']
    x1 = dm['Coordinates'][:,0]#/1e3/(1.0+zl) # Mpc/h, comoving
    x2 = dm['Coordinates'][:,1]#/1e3/(1.0+zl) # Mpc/h, comoving
    x3 = dm['Coordinates'][:,2]#/1e3/(1.0+zl) # Mpc/h, comoving

    ptnl_p = dm['Potential']
    idx_ptnl_p_min = ptnl_p == ptnl_p.min()

    xcp1 = x1[idx_ptnl_p_min]
    xcp2 = x2[idx_ptnl_p_min]
    xcp3 = x3[idx_ptnl_p_min]

    #sdens = pcs.call_sph_sdens_arrays(x1,x2,x3,xc1,xc2,xc3,bsz,nnn,npp,pmass)*1e6
    sdens = read_sdens_from_illustris(dm,x1-xcp1,x2-xcp2,x3-xcp3,bsz,nnn)*(1.0+zl)**2.0
    #kappa = sdens/mm.sigma_crit(zl,zs0)

    #pl.figure()
    #pl.contourf(sdens)
    #pl.colorbar()
    ##print np.max(kappa)
    #----------------------------------------------------------------------
    bsz = bsz/1e3/mm.Da(zl)*mm.apr #arcsec
    dsx = bsz/nnn #arcsec

    print bsz

    xi1,xi2  = make_r_coor(nnn,dsx)
    ai1,ai2 = pcs.call_cal_alphas0(sdens,nnn,bsz,zl,zs0)
    lensed_img = pcs.call_ray_tracing_single(ai1,ai2,nnn,bsz,zl)
    lenses_img = slf.loadin_fits(HaloID,zl,bsz,nnn)

#----------------------------------------------------------------------
    #zs = 2.0
    #dsi = 0.03
    #ysc1 = 0.0
    #ysc2 = 0.0
    #ai1,ai2 = pcs.call_cal_alphas0(sdens,nnn,bsz,zl,zs)
    #yi1 = xi1-ai1
    #yi2 = xi2-ai2

    #g_source = pyfits.getdata("./dry_run_pipeline/439.0_149.482739_1.889989_processed.fits")
    #g_source = np.array(g_source,dtype="<d")
    #g_source[g_source<=0.0001] = 1e-6
    #lensed_img = lv4.call_ray_tracing(g_source,yi1,yi2,ysc1,ysc2,dsi)
#----------------------------------------------------------------------
    pl.figure()
    pl.contourf(lensed_img)
    pl.colorbar()

    final_img = lensed_img+lenses_img/100.0

    pyfits.writeto("test.fits",final_img,clobber=True)

    pl.figure()
    pl.contourf(final_img)
    pl.colorbar()
    return 0
if __name__ == '__main__':
    run_main()
    pl.show()
