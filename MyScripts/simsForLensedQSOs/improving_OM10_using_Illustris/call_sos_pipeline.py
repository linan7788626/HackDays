#/usr/bin/env python
import sys
del sys.path[2]
import numpy as np
#import astropy.io.fits as pyfits
import pylab as pl
import python_call_c_so as pcs
import triangle_root_finding as trf
import illustris_python as il
import mycosmology as mm
import scipy.interpolate as sci
import matplotlib as mpl
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
def point_ray_tracing(xi1,xi2,ai1,ai2,spar,bsz,nnn):
    g_limage = xi1*0.0
    dsx = bsz/nnn

    xroot1,xroot2,nroots = trf.roots_zeros(xi1,xi2,ai1,ai2,spar[0],spar[1])

    idr1 = ((np.array(xroot1)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    idr2 = ((np.array(xroot2)+bsz/2.0-dsx/2.0)/dsx).astype('int')
    g_limage[idr1,idr2] = spar[2]

    return g_limage

def calculate_new_ys(ximg1,ximg2,alpha1,alpha2,dsi):

    aimg1 = pcs.call_inverse_cic(alpha1,ximg1,ximg2,dsi)
    aimg2 = pcs.call_inverse_cic(alpha2,ximg1,ximg2,dsi)

    yimg1 = ximg1-aimg1
    yimg2 = ximg2-aimg2

    ysf1 = np.mean(yimg1)
    ysf2 = np.mean(yimg2)

    return ysf1,ysf2
#------------------------------------------------------------------------------
def run_main():
    #----------------------------------------------------------------------
    SnapID = 135
    HaloID = 50

    zl = 0.5
    zs = 2.0
    nnn=512

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

    bsz=50.0 # kpc/h
    sdens = read_sdens_from_illustris(dm,x1-xcp1,x2-xcp2,x3-xcp3,bsz,nnn)
    #kappa = sdens/mm.sigma_crit(zl,zs)
    #----------------------------------------------------------------------
    bsz = bsz/1e3/mm.Da(zl)*mm.apr #arcsec
    dsx = bsz/nnn #arcsec

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

    ldays = np.zeros((nroots,len(time_days)))
    #for i in xrange(nroots):
        #ldays[i,:] = time_days+td[i]
    ldays[0,:] = time_days+td[0]
    ldays[1,:] = time_days+td[1]
    ldays[2,:] = time_days+td[2]

    mags = mags.astype("double")

    lmags = np.zeros((nroots,len(time_days)))
    #for i in xrange(nroots):
        #lmags[i,:] = mags-2.5*np.log10(np.abs(muii[i]))
    lmags[0,:] = mags-2.5*np.log10(np.abs(muii[0]))
    lmags[1,:] = mags-2.5*np.log10(np.abs(muii[1]))
    lmags[2,:] = mags-2.5*np.log10(np.abs(muii[2]))

    cmap = mpl.cm.jet_r

#------------------------------------------------------------------------------
#linestyle=linestyle, color=color
    pl.figure(figsize=(8,8))
    pl.xlabel("Days")
    pl.ylabel("Mags")
    pl.ylim(27,21)
    #for i in xrange(nroots):
        #pl.plot(ldays[i,:]-49500+820,lmags[i,:],linestyle='-',color=cmap(i/float(nroots)))
    pl.plot(ldays[0,:]-49500+820,lmags[0,:],linestyle='-',color="r")
    pl.plot(ldays[1,:]-49500+820,lmags[1,:],linestyle='-',color="g")
    pl.plot(ldays[2,:]-49500+820,lmags[2,:],linestyle='-',color="b")

    pl.figure(figsize=(8,8))
    pl.xlim(-5.0,5.0)
    pl.ylim(-5.0,5.0)
    pl.xlabel("arcsec")
    pl.ylabel("arcsec")
    #for i in xrange(nroots):
        #pl.plot(xroot1[i],xroot2[i],marker='o',color=cmap(i/float(nroots)))
#markersize=10
    pl.plot(xroot1[0],xroot2[0],marker='o',color="r",markersize=10)
    pl.plot(xroot1[1],xroot2[1],marker='o',color="g",markersize=10)
    pl.plot(xroot1[2],xroot2[2],marker='o',color="b",markersize=10)

    pl.contour(xi1,xi2,mua,colors=('r',))
    pl.contour(yi1,yi2,mua,colors=('g',))
    pl.plot(ys1,ys2,'ks')

    print td
if __name__ == '__main__':
    run_main()
    pl.show()
