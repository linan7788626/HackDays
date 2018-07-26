import numpy as np
#from pylab import *
#import scipy.signal as ss
import mycosmology as mm
import libv4_cv as lv4
import pyfits
import sys

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
#--------------------------------------------------------------------
zl0 = 0.1
zs = 5.0
sigmav0 = 320			#km/s
q0 = 0.0
rc0 = 0.0

boxsize = 6.0 # (arcsec)
nnn = 512
dsx = boxsize/nnn
dsi = dsx
npl = 3
zln = np.linspace(0.0,zs,npl+2)[1:-1]

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
sigmav[1] = sigmav0
#sigmav[-1] = sigmav0

xi1[0] = xx1
xi2[0] = xx2
if (npl==1):
	ag1[0],ag2[0] = alpha_nie(sigmav[0],xx1,xx2,q0,rc0,zln[0],zs)
	ai1[0] = itpn_on_kth_plane(ag1[0],xi1[0],xi2[0])
	ai2[0] = itpn_on_kth_plane(ag2[0],xi1[0],xi2[0])

	af1 = ai1[0]
	af2 = ai2[0]

if (npl==2):
	ag1[0],ag2[0] = alpha_nie(sigmav[0],xx1,xx2,q0,rc0,zln[0],zs)
	ai1[0] = itpn_on_kth_plane(ag1[0],xi1[0],xi2[0])
	ai2[0] = itpn_on_kth_plane(ag2[0],xi1[0],xi2[0])

	B01 = Bij(zln[0],zln[1])
	xi1[1] = xi1[0]-ai1[0]*B01
	xi2[1] = xi2[0]-ai2[0]*B01
	ag1[1],ag2[1] = alpha_nie(sigmav[1],xx1,xx2,q0,rc0,zln[1],zs)
	ai1[1] = itpn_on_kth_plane(ag1[1],xi1[1],xi2[1])
	ai2[1] = itpn_on_kth_plane(ag2[1],xi1[1],xi2[1])

	af1 = ai1[0]+ai1[1]
	af2 = ai2[0]+ai2[1]

if (npl>2):
	ag1[0],ag2[0] = alpha_nie(sigmav[0],xx1,xx2,q0,rc0,zln[0],zs)
	ai1[0] = itpn_on_kth_plane(ag1[0],xi1[0],xi2[0])
	ai2[0] = itpn_on_kth_plane(ag2[0],xi1[0],xi2[0])

	B01 = Bij(zln[0],zln[1])
	xi1[1] = xi1[0]-ai1[0]*B01
	xi2[1] = xi2[0]-ai2[0]*B01
	ag1[1],ag2[1] = alpha_nie(sigmav[1],xx1,xx2,q0,rc0,zln[1],zs)
	ai1[1] = itpn_on_kth_plane(ag1[1],xi1[1],xi2[1])
	ai2[1] = itpn_on_kth_plane(ag2[1],xi1[1],xi2[1])

	af1 = ai1[0]+ai1[1]
	af2 = ai2[0]+ai2[1]
	for j in xrange(2,npl):
		sum_alpha1 = yi1*0.0
		sum_alpha2 = yi2*0.0
		for i in xrange(j):
			Btmp1 = Bij(zln[i],zln[j])
			Btmp2 = Bij(zln[i],zln[j-1])
			sum_alpha1 += (Btmp1-Btmp2)*ai1[i]
			sum_alpha2 += (Btmp1-Btmp2)*ai2[i]

		xi1[j] = xi1[j-1]-sum_alpha1
		xi2[j] = xi2[j-1]-sum_alpha2
		ag1[j],ag2[j] = alpha_nie(sigmav[j],xx1,xx2,q0,rc0,zln[j],zs)
		ai1[j] = itpn_on_kth_plane(ag1[j],xi1[j],xi2[j])
		ai2[j] = itpn_on_kth_plane(ag2[j],xi1[j],xi2[j])

		af1 += ai1[j]
		af2 += ai2[j]

yi1 =  xi1[0]-af1
yi2 =  xi2[0]-af2

#----------------------------------------------------------------------
g_amp = 1.0   	# peak brightness value
g_sig = 0.01  	# Gaussian "sigma" (i.e., size)
g_xcen = 0.0  	# x position of center (also try (0.0,0.14)
g_ycen = 0.0  	# y position of center
g_axrat = 1.0 	# minor-to-major axis ratio
g_pa = 0.0    	# major-axis position angle (degrees) c.c.w. from x axis
#----------------------------------------------------------------------
#g_source = 0.0*xx1
#gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
#g_source = gauss_2d(xx1,xx2,gpar)
#g_lensimage = 0.0*yi1
#g_lensimage = gauss_2d(yi1,yi2,gpar)

#g_source = pyfits.getdata("./compound_R_1_0_S_1_1_u.fits")
input_fits = sys.argv[1]
g_source = pyfits.getdata(input_fits)
g_source = np.array(g_source,dtype="<d")
g_source_pin = lv4.call_ray_tracing(g_source,xx1,xx2,g_xcen,g_ycen,dsi)
g_lensimage = lv4.call_ray_tracing(g_source,yi1,yi2,g_xcen,g_ycen,dsi)


filename = "output.fits"
pyfits.writeto(filename,g_lensimage,clobber=True)
##--------------------------lens images contour------------------------
##levels = [0.15,0.30,0.45,0.60,0.75,0.9,1.05]
#figure(num=None,figsize=(10,5),dpi=80, facecolor='w', edgecolor='k')


#a = axes([0.05,0.1,0.4,0.8])
#a.set_xlim(-boxsize/4.0,boxsize/4.0)
#a.set_ylim(-boxsize/4.0,boxsize/4.0)
#a.contourf(xx1,xx2,g_source_pin)
##a.contour(yi1,yi2,mua,0,colors=('g'),linewidths = 2.0)


#b = axes([0.55,0.1,0.4,0.8])
#b.set_xlim(-boxsize/4.0,boxsize/4.0)
#b.set_ylim(-boxsize/4.0,boxsize/4.0)
#b.contourf(xx1,xx2,g_lensimage)
##b.contour(xi1,xi2,muap,colors=('k'),linewidths = 2.0)
##savefig('output.pdf')
#show()
#------------------------------------------------------------------------------
