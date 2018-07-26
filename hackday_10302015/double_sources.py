import numpy as np
#from pylab import *
#import scipy.signal as ss
import mycosmology as mm
import libv4_cv as lv4
import pyfits
import sys
import pylab as pl

import scipy.optimize as sco

def parabola_1d(x,xb,xc,a):
    res = a*(x-xb)*(2.0*xc-(x-xb))
    res[res<=0] = 0.0
    return res

def isNewImage(x1,x2,xil1,xil2):
    rdist = np.hypot((x1-xil1),(x2-xil2))
    lidx = len(rdist[rdist < 1e-7])
    return lidx

def root_finding(x_guess,y10,y20):
    def simple_lensing_equation(x):
        y = x*0.0
        #y[0],y[1] = nie_all(x[0],x[1],xlc1,xlc2,re0,rc0,ql0,phi0,y10,y20)[9:]
        y[0],y[1] = leq_nie_zeros(320.0,x[0],x[1],0.05,0.00,0.5,1.0,y10,y20);
        return [y[0],y[1]]

    sol = sco.root(simple_lensing_equation,[x_guess[0],x_guess[1]],method='krylov')
    #,options={'xtol':1e-9})
    # 'anderson','hybr','lm','broyden1','broyden2','anderson','linearmixing','diagbroyden','excitingmixing','krylov','df-sane'
    return sol.x

#--------------------------------------------------------------------
def re_sv(sv,z1,z2):
    #Dl = mm.Da(z1)
    Ds = mm.Da(z2)
    Dds = mm.Da2(z1,z2)

    res = 4.0*np.pi*(sv**2.0/mm.vc**2.0)*Dds/Ds*mm.apr
    return res

#--------------------------------------------------------------------
def img_cic(ysc2_in,ysc1_in,dsi_in,img_in,x2_in,x1_in):

    img_out=img_in*0.0

    nx = np.shape(img_in)[0]
    ny = np.shape(img_in)[1]

    xc1=ysc1_in
    xc2=ysc2_in

    xb1 = (x1_in-xc1)/dsx+nx/2.0-0.5
    xb2 = (x2_in-xc2)/dsx+ny/2.0-0.5

    #i1 = xb1.astype(int)
    #j1 = xb2.astype(int)
    #wx = 1.-(xb1-i1)
    #wy = 1.-(xb2-j1)
    #iarr = np.array([i1,i1,i1+1,i1+1])
    #jarr = np.array([j1,j1+1,j1,j1+1])
    #warr = np.array([wx*wy,wx*(1.0-wy),(1.0-wx)*wy,(1.0-wx)*(1.0-wy)])
    #img_out_tmp = img_in[iarr,jarr]*warr
    #img_out = np.sum(img_out_tmp,axis=0)
    #img_out = img_in[i1,j1]
    i1 = xb1.astype(int)
    wx = 1.-(xb1-i1)

    j1 = xb2.astype(int)
    wy = 1.-(xb2-j1)

    ww1 = wx*wy
    ww2 = wx*(1.0-wy)
    ww3 = (1.0-wx)*wy
    ww4 = (1.0-wx)*(1.0-wy)

    idxi1 = i1>nx-2
    idxi2 = i1<0

    idxj1 = j1>ny-2
    idxj2 = j1<0

    idx = idxi1|idxi2|idxj1|idxj2
    ww1[idx] = 0.0
    ww2[idx] = 0.0
    ww3[idx] = 0.0
    ww4[idx] = 0.0

    i1[idxi2] = 0
    i1[idxi1] = nx-2
    j1[idxj2] = 0
    j1[idxj1] = ny-2

    iarr = np.array([i1,i1,i1+1,i1+1])
    jarr = np.array([j1,j1+1,j1,j1+1])
    warr = np.array([ww1,ww2,ww3,ww4])
    img_out_tmp = img_in[iarr,jarr]*warr

    img_out = np.sum(img_out_tmp,axis=0)
    #print np.max(img_in),np.max(img_out)

    return img_out
#--------------------------------------------------------------------
def itpn_on_kth_plane(img_in,x2_in,x1_in):

    img_out=img_in*0.0

    nx = np.shape(img_in)[0]
    ny = np.shape(img_in)[1]

    xc1=0.0
    xc2=0.0

    dsx_in = abs(x1_in[1,0]-x1_in[0,0])
    dsy_in = abs(x2_in[0,1]-x2_in[0,0])

    xb1 = (x1_in-xc1)/dsx_in+nx/2.0-0.5
    xb2 = (x2_in-xc2)/dsy_in+ny/2.0-0.5

    i1 = xb1.astype(int)
    wx = 1.-(xb1-i1)

    j1 = xb2.astype(int)
    wy = 1.-(xb2-j1)

    ww1 = wx*wy
    ww2 = wx*(1.0-wy)
    ww3 = (1.0-wx)*wy
    ww4 = (1.0-wx)*(1.0-wy)

    idxi1 = i1>nx-2
    idxi2 = i1<0

    idxj1 = j1>ny-2
    idxj2 = j1<0

    idx = idxi1|idxi2|idxj1|idxj2
    ww1[idx] = 0.0
    ww2[idx] = 0.0
    ww3[idx] = 0.0
    ww4[idx] = 0.0

    i1[idxi2] = 0
    i1[idxi1] = nx-2
    j1[idxj2] = 0
    j1[idxj1] = ny-2

    iarr = np.array([i1,i1,i1+1,i1+1])
    jarr = np.array([j1,j1+1,j1,j1+1])
    warr = np.array([ww1,ww2,ww3,ww4])
    img_out_tmp = img_in[iarr,jarr]*warr

    img_out = np.sum(img_out_tmp,axis=0)
    #print np.max(img_in),np.max(img_out)

    return img_out
#--------------------------------------------------------------------
def kappa_nie(sigmav,x1,x2,qe,xcore,z1,z2):
    scr = mm.sigma_crit(z1,z2)
    x = np.sqrt(xcore*xcore+x1*x1*(1.0-qe)+x2*x2*(1.0+qe))
    ac = sigmav**2.0*mm.apr/(mm.G*mm.Da(z1)*scr)
    res = ac/(2.0*x)
    return res
def phi_nie(sigmav,x1,x2,qe,xcore,z1,z2):
    scr = mm.sigma_crit(z1,z2)
    x = np.sqrt(xcore*xcore+x1*x1*(1.0-qe)+x2*x2*(1.0+qe))
    #res2 = Dds/Ds*4.0*np.pi*sigmav**2.0/vc**2.0*mm.apr
    ac = sigmav**2.0*mm.apr/(mm.G*mm.Da(z1)*scr)
    res = ac*x
    return res

def alpha_nie(sigmav,x1,x2,qe,xcore,z1,z2):
    scr = mm.sigma_crit(z1,z2+1e-7)
    x = np.sqrt(xcore*xcore+x1*x1*(1.0-qe)+x2*x2*(1.0+qe))
    ac = sigmav**2.0*mm.apr/(mm.G*mm.Da(z1)*scr)
    res1 = ac*(1-qe)*x1/x
    res2 = ac*(1+qe)*x2/x
    return res1,res2

def leq_nie_zeros(sigmav,x1,x2,qe,xcore,z1,z2,y10,y20):
    scr = mm.sigma_crit(z1,z2+1e-7)
    x = np.sqrt(xcore*xcore+x1*x1*(1.0-qe)+x2*x2*(1.0+qe))
    ac = sigmav**2.0*mm.apr/(mm.G*mm.Da(z1)*scr)
    res1 = ac*(1-qe)*x1/x
    res2 = ac*(1+qe)*x2/x
    return x1-res1-y10,x2-res2-y20

def green_iso(x,y):
    epsilon= 0.000001
    r = np.sqrt(x**2.0+y**2.0+epsilon**2.0)
    res = 1.0/np.pi*np.log(r)
    return res
#--------------------------------------------------------------------
def lensing_signals_phi(phi_in,dbx):
    phi2,phi1 = np.gradient(phi_in,dbx)
    phi12,phi11 = np.gradient(phi1,dbx)
    phi22,phi21 = np.gradient(phi2,dbx)

    alpha1 = phi1
    alpha2 = phi2
    kappa  = (phi11+phi22)*0.5
    shear1 = (phi11-phi22)*0.5
    shear2 = (phi12+phi21)*0.5

    mu = 1.0/(1.0-phi11-phi22+phi11*phi22-phi12*phi21)

    return alpha1,alpha2,kappa,shear1,shear2,mu
#--------------------------------------------------------------------
def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x-xcen)*np.cos(phirad)+(y-ycen)*np.sin(phirad)
    ynew = (y-ycen)*np.cos(phirad)-(x-xcen)*np.sin(phirad)
    return (xnew,ynew)
def gauss_2d(x, y, par):
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    res0 = ((xnew**2)*par[4]+(ynew**2)/par[4])/np.abs(par[1])**2
    res = par[0]*np.exp(-0.5*res0)
    return res
def Bij(z1,z2):
    res = mm.Da2(z1,z2)*mm.Da(zs)/(mm.Da(z2)*mm.Da2(z1,zs))
    return res

def td_nie(xi1,xi2,q0,rc0,z1,z2):

    dsx_in = xi1[1,1]-xi1[0,0]

    ai1,ai2 = alpha_nie(sigmav[0],xx1,xx2,q0,rc0,z1,z2)
    phi =  phi_nie(sigmav[0],xx1,xx2,q0,rc0,z1,z2)
    #Kc = 1.0
    Kc = (1.0+zl)/mm.vc*(mm.Da(z1)*mm.Da(z2)/mm.Da2(z1,z2))*mm.Mpc/(1e3*mm.dy)/mm.apr/mm.apr
    td = Kc*(0.5*((ai1)**2.0+(ai2)**2.0)-phi)

    print Kc
    print np.max(td),
    print np.min(td)
    print np.max(td)-np.min(td)

    phi12,phi11 = np.gradient(ai1,dsx_in)
    phi22,phi21 = np.gradient(ai2,dsx_in)

    mu = 1.0/(1.0-phi11-phi22+phi11*phi22-phi12*phi21)

    return td,mu
#--------------------------------------------------------------------
zl0 = 0.1
zl = 0.5
zs0 = 10.0
zs1 = 1.0
zs2 = 5.0
zs = 2.0
sigmav0 = 320           #km/s
q0 = 0.05
rc0 = 0.0

boxsize = 10.0 # (arcsec)
nnn = 128
dsx = boxsize/nnn
dsi = 1.0*dsx
npl = 3
zln = np.linspace(0.0,zs,npl+2)[1:-1]
zln[0] = zl

xx01 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
xx02 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
xx1,xx2 = np.meshgrid(xx01,xx02)

xi1 = np.zeros((npl,nnn,nnn))
xi2 = np.zeros((npl,nnn,nnn))
ag1 = np.zeros((npl,nnn,nnn))
ag2 = np.zeros((npl,nnn,nnn))
ai1 = np.zeros((npl,nnn,nnn))
ai2 = np.zeros((npl,nnn,nnn))
yi1 = np.zeros((nnn,nnn))
yi2 = np.zeros((nnn,nnn))

sigmav = np.zeros((npl))
sigmav[0] = sigmav0
af11,af12 = alpha_nie(sigmav[0],xx1,xx2,q0,rc0,zl,zs1)

yi11 =  xx1-af11#*(mm.Da2(zl,zs1)/mm.Da(zs1))*(mm.Da(zs0)/mm.Da2(zl,zs0))
yi12 =  xx2-af12#*(mm.Da2(zl,zs1)/mm.Da(zs1))*(mm.Da(zs0)/mm.Da2(zl,zs0))

af21,af22 = alpha_nie(sigmav[0],xx1,xx2,q0,rc0,zl,zs2)

yi21 =  xx1-af21#*(mm.Da2(zl,zs2)/mm.Da(zs2))*(mm.Da(zs0)/mm.Da2(zl,zs0))
yi22 =  xx2-af22#*(mm.Da2(zl,zs2)/mm.Da(zs2))*(mm.Da(zs0)/mm.Da2(zl,zs0))

#----------------------------------------------------------------------
g_amp = 1.0     # peak brightness value
g_sig = 0.01    # Gaussian "sigma" (i.e., size)
g_xcen = 0.0    # x position of center (also try (0.0,0.14)
g_ycen = 0.0    # y position of center
g_axrat = 1.0   # minor-to-major axis ratio
g_pa = 0.0      # major-axis position angle (degrees) c.c.w. from x axis
#----------------------------------------------------------------------

def output_lensed_images(input_fits,yy1,yy2,ys1,ys2,dsi):
    g_source = pyfits.getdata(input_fits)
    g_source = np.array(g_source,dtype="<d")
    #g_source_pin = lv4.call_ray_tracing(g_source,xx1,xx2,g_xcen,g_ycen,dsi)
    g_lensimage = lv4.call_ray_tracing(g_source,yy1,yy2,ys1,ys2,dsi)
    return g_lensimage


ysn1 = 0.01
ysn2 = 0.02
re_0 = re_sv(sigmav0,zl,zs2)

xr1 = np.zeros((5))
xr2 = np.zeros((5))
xg1 = np.array([ysn1-re_0,ysn1,ysn1,ysn1+re_0])
xg2 = np.array([ysn2,ysn2-re_0,ysn2+re_0,ysn2])
ncount = 0

g_lsn = np.zeros((nnn,nnn))

for i in xrange(len(xg1)):
    xrt1,xrt2 = root_finding([xg1[i],xg2[i]],ysn1,ysn2)
    if isNewImage(xrt1,xrt2,xr1,xr2) <= 0:
        xr1[ncount]=xrt1
        xr2[ncount]=xrt2
        ncount = ncount + 1


xr1_idx = (xr1+boxsize/2.0-dsx/2.0)/dsx

xr1_idx = xr1_idx.astype("int")
xr2_idx = (xr2+boxsize/2.0-dsx/2.0)/dsx
xr2_idx = xr2_idx.astype("int")

g_lsn[xr1_idx,xr2_idx] = 0.01

ys11 = 0.0
ys12 = 0.0
ys21 = 0.0
ys22 = 0.0

input1 = sys.argv[1]
input2 = sys.argv[2]
image1 = output_lensed_images(input1,yi11,yi12,ys11,ys12,dsi)
image2 = output_lensed_images(input2,yi21,yi22,ys21,ys22,dsi)

td_sn,mu_sn = td_nie(xx1,xx2,0.05,0.0,zl,zs1)

#print np.max(td_sn),np.max(mu_sn)

#print np.max(image1),np.max(image2)

ic = 10.0

#levels = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6]
levels = [0.0,0.2,1.2,2.4,3.6]
time_days,mags = np.loadtxt("./SN_opsim.csv",dtype='string', delimiter=',', usecols=(1, 4), unpack=True)
time_days = np.double(time_days).astype("int")
mags = mags.astype("double")

#print len(time_days)
#time_days_sub = time_days[20::5]
#mags_sub = mags[20::5]
#nsub = len(time_days)/npics
#print len(time_days_sub)

print time_days

mags_gals = 17.0+9.0
nnan = len(mags[np.isnan(mags)])
mags_sns = mags[nnan:]
rat = 10.0**(mags_gals-mags_sns)
npics = len(mags_sns)
#print mags_sns

#print mags_sns
#print rat

for i in xrange(npics):
    sktd = td_sn/td_sn.max()*20.0
    idx = time_days[i]+sktd-time_days[0]
    idx = idx.astype("int")
    ratio_map = rat[idx]

    final_image = image1+image2+g_lsn*np.abs(mu_sn)*ratio_map
    #print np.max(g_lsn*np.abs(mu_sn)*ratio_map)
    #print np.max(final_image)

    #pl.imshow(final_image,cmap = pl.get_cmap(gray))
    #pl.figure(figsize=(10,10))
    #pl.contourf(final_image,levels)
    #pl.savefig("./output_pngs/"+'{:0>10}'.format(str(npics-i))+".png")

    filename = "./output_fits/"+'{:0>10}'.format(str(i))+"_output_double.fits"
    #pyfits.writeto(filename,image1+image2+g_lsn*np.abs(mu_sn),clobber=True)
    pyfits.writeto(filename,final_image,clobber=True)

##filename = "output_double.fits"
##pyfits.writeto(filename,image1+image2+g_lsn*np.abs(mu_sn),clobber=True)
###--------------------------lens images contour------------------------
###levels = [0.15,0.30,0.45,0.60,0.75,0.9,1.05]
##figure(num=None,figsize=(10,5),dpi=80, facecolor='w', edgecolor='k')


##a = axes([0.05,0.1,0.4,0.8])
##a.set_xlim(-boxsize/4.0,boxsize/4.0)
##a.set_ylim(-boxsize/4.0,boxsize/4.0)
##a.contourf(xx1,xx2,g_source_pin)
###a.contour(yi1,yi2,mua,0,colors=('g'),linewidths = 2.0)


##b = axes([0.55,0.1,0.4,0.8])
##b.set_xlim(-boxsize/4.0,boxsize/4.0)
##b.set_ylim(-boxsize/4.0,boxsize/4.0)
##b.contourf(xx1,xx2,g_lensimage)
###b.contour(xi1,xi2,muap,colors=('k'),linewidths = 2.0)
###savefig('output.pdf')
#pl.show()
##------------------------------------------------------------------------------
