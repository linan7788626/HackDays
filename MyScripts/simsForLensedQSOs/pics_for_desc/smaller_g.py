#/usr/bin/env python
import sys
del sys.path[2]
import numpy as np
import astropy.io.fits as pyfits
import pylab as pl
import python_call_c_so as pcs
import triangle_root_finding as trf
import illustris_python as il
import mycosmology as mm
import scipy.interpolate as sci
import matplotlib as mpl
# import subID_loadin_fits as slf
import libv4_cv as lv4
basePath = '../Illustris-3'
##---------------------------------------------------------------------------------
#def read_sdens_from_illustris(xn1,xn2,xn3,bsz,nnn):
    #dsx = bsz/nnn

    #points = np.vstack([xn1,xn2,xn3]).T
    #values = dm['SubfindDensity']*1e10
    #grid_x, grid_y, grid_z = np.mgrid[0:1:(nnn*1j),0:1:(nnn*1j),0:1:(nnn*1j)]*bsz-bsz/2.0+dsx/2.0
    #xi1,xi2 = np.mgrid[0:1:(nnn*1j),0:1:(nnn*1j)]*bsz-bsz/2.0+dsx/2.0

    #dens_3d = sci.griddata(points,values,(grid_x,grid_y,grid_z),method='linear',fill_value=0)
    #sdens = np.sum(dens_3d,axis=2)*dsx*1e6
    #return sdens

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
def stacking_points_and_galaxies(xroot1,xroot2,lroots,lensed_imgs,bsz,nnn):
    res = lensed_imgs*0.0
    dsx = bsz/nnn

    idr1 = ((np.array(xroot1)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    idr2 = ((np.array(xroot2)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    res[idr1,idr2] = lroots

    res = res + lensed_imgs

    return res
#------------------------------------------------------------------------------
def run_main():
    #----------------------------------------------------------------------
    SnapID = 135
    HaloID = 50

    bsz=200.0 # kpc/h
    zl = 0.5
    zs = 2.0
    #zs0= 10.0
    nnn=512
    dsi = 0.03

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

    sdens = read_sdens_from_illustris(dm,x1-xcp1,x2-xcp2,x3-xcp3,bsz,nnn)
    #kappa = sdens/mm.sigma_crit(zl,zs)
    #----------------------------------------------------------------------
    bsz = bsz/1e3/mm.Da(zl)*mm.apr #arcsec
    dsx = bsz/nnn #arcsec

    print bsz
    print dsx

    xi1,xi2  = make_r_coor(nnn,dsx)
    ai1,ai2 = pcs.call_cal_alphas0(sdens,nnn,bsz,zl,zs)
    yi1 = xi1 - ai1
    yi2 = xi2 - ai2
    phi = pcs.call_cal_phi(sdens,nnn,bsz,zl,zs)
    mua = pcs.call_alphas_to_mu(ai1,ai2,nnn,dsx)
    #----------------------------------------------------------------------
    ys1 = 0.0
    ys2 = 0.0
    xroot1,xroot2,nroots = trf.roots_zeros(xi1,xi2,ai1,ai2,ys1,ys2)

    #pl.figure()
    #pl.plot(ys1,ys2,'ko')
    #pl.plot(xroot1,xroot2,'rs')
    #pl.contour(xi1,xi2,mua)
    #pl.contour(yi1,yi2,mua)
    #pl.colorbar()

    muii = pcs.call_inverse_cic_single(mua,xroot1,xroot2,dsx)
    phii = pcs.call_inverse_cic_single(phi,xroot1,xroot2,dsx)

    Kc = (1.0+zl)/mm.vc*(mm.Da(zl)*mm.Da(zs)/mm.Da2(zl,zs))*mm.Mpc/(1e3*mm.dy)/mm.apr/mm.apr
    td = Kc*(0.5*((xroot1-ys1)**2.0+(xroot2-ys2)**2.0)-phii)

    time_days,mags = np.loadtxt("./SN_opsim.csv",dtype='string', delimiter=',', usecols=(1, 4), unpack=True)
    time_days = time_days.astype("double")
#------------------------------------------------------------------------------
    ldays = np.zeros((nroots,len(time_days)))
    for i in xrange(nroots):
        ldays[i,:] = time_days+td[i]

    mags = mags.astype("double")
    lmags = np.zeros((nroots,len(time_days)))

    for i in xrange(nroots):
        lmags[i,:] = mags-2.5*np.log10(np.abs(muii[i]))

    #cmap = mpl.cm.jet_r
    cmap = mpl.cm.gist_earth

    pl.figure(figsize=(8,8))
    pl.xlabel("Days")
    pl.ylabel("Mags")
    pl.ylim(30,22)
    for i in xrange(nroots):
        pl.plot(ldays[i,:]-49500+31010-50,lmags[i,:],linestyle='-',color=cmap(i/float(nroots)),markersize=10)
    #pl.plot(ldays[0,:]-49500+31010-50,lmags[0,:],linestyle='-',color="k")
    #pl.plot(ldays[1,:]-49500+31010-50,lmags[1,:],linestyle='-',color="b")
    #pl.plot(ldays[2,:]-49500+31010-50,lmags[2,:],linestyle='-',color="g")
    #pl.plot(ldays[3,:]-49500+31010-50,lmags[3,:],linestyle='-',color="m")
    #pl.plot(ldays[4,:]-49500+31010-50,lmags[4,:],linestyle='-',color="r")

    pl.figure(figsize=(8,8))
    #pl.xlim(-bsz/2.0,bsz/2.0)
    #pl.ylim(-bsz/2.0,bsz/2.0)
    pl.xlim(-10.0,10.0)
    pl.ylim(-10.0,10.0)
    pl.xlabel("arcsec")
    pl.ylabel("arcsec")
    for i in xrange(nroots):
        pl.plot(xroot1[i],xroot2[i],marker='o',color=cmap(i/float(nroots)),markersize=10)
    #pl.plot(xroot1[0],xroot2[0],marker='o',color="k",markersize=15) #
    #pl.plot(xroot1[1],xroot2[1],marker='o',color="b",markersize=15) #
    #pl.plot(xroot1[2],xroot2[2],marker='o',color="g",markersize=15) #
    #pl.plot(xroot1[3],xroot2[3],marker='o',color="m",markersize=15) #
    #pl.plot(xroot1[4],xroot2[4],marker='o',color="r",markersize=15) #

    pl.contour(xi1,xi2,mua,colors=('r',))
    pl.contour(yi1,yi2,mua,colors=('g',))
    pl.plot(ys1,ys2,'cs')
    print td
#------------------------------------------------------------------------------
    ysc1 = 0.0
    ysc2 = 0.0

    g_source = pyfits.getdata("./dry_run_pipeline/439.0_149.482739_1.889989_processed.fits")
    g_source = np.array(g_source,dtype="<d")
    pl.figure()
    pl.contourf(g_source)
    pl.colorbar()
    std_srcs = np.std(g_source)
    print std_srcs,std_srcs*4.0
    gthhd = 6.0*std_srcs
    g_source = g_source - gthhd
    g_source[g_source<=0.0] = 0.0
    lensed_img = lv4.call_ray_tracing(g_source,yi1,yi2,ysc1,ysc2,dsi)
    #lensed_img = pcs.call_ray_tracing_single(ai1,ai2,nnn,bsz,zl)
    # lenses_img = slf.loadin_fits(HaloID,zl,bsz,nnn)

    lroots =10.0**((lmags-25.523)/(-2.5))

    final_img = lensed_img#+lenses_img/100.0
    res = stacking_points_and_galaxies(xroot1,xroot2,lroots[:,70],final_img,bsz,nnn)
    pl.figure()
    pl.contourf(xi1,xi2,res)
    pl.colorbar()
    pyfits.writeto("all.fits",res,clobber=True)
    return res

if __name__ == '__main__':
    run_main()
    pl.show()
