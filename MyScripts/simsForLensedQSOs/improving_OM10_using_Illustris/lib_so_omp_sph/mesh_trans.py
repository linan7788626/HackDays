import numpy as np
import pylab as pl
import scipy.ndimage.filters as snf
#--------------------------------------------------------------------
zl = 0.1
zs = 1.0
sigmav = 220			#km/s
q0 = 0.2
rc0 = 0.1
##--------------------------------------------------------------------
#def re_sv(sv,z1,z2):
#	Dl = Da(z1)
#	Ds = Da(z2)
#	Dds = Da2(z1,z2)
#
#	res = 4.0*np.pi*(sv**2.0/vc**2.0)*Dds/Ds*apr
#	return res
#--------------------------------------------------------------------
def lens_equation_sie(x1,x2,lpar):
    xc1 = lpar[0]   #x coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]   #y coordinate of the center of lens (in units of Einstein radius).
    q   = lpar[2]   #Ellipticity of lens.
    rc  = lpar[3]   #Core size of lens (in units of Einstein radius).
    re  = lpar[4]   #Einstein radius of lens.
    pha = lpar[5]   #Orintation of lens.

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

    y11 = 1.0-rea11
    y22 = 1.0-rea22
    y12 = 0.0-rea12
    y21 = 0.0-rea21

    jacobian = y11*y22-y12*y21
    mu = 1.0/jacobian

    res1 = (a1*cosa-a2*sina)*re
    res2 = (a2*cosa+a1*sina)*re
    return res1,res2,mu
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

def pixel_trans(x1,x2,xx1,xx2,matrix,ntmp):

    ntmp = ntmp*1
    dr = 1000
    nnx,nny = np.shape(matrix)
    dsx = xx1[1,1]-xx1[0,0]
    bsz = nnx*dsx

    buf = np.zeros((2.0*ntmp,2.0*ntmp))
    buf_x1 = np.zeros((2.0*ntmp,2.0*ntmp))
    buf_x2 = np.zeros((2.0*ntmp,2.0*ntmp))

    kk = 0

    while (dr >dsx) :
        i = int((x2+bsz*0.5)/dsx)
        j = int((x1+bsz*0.5)/dsx)

        idx1 = i-ntmp;sidx1 = 0
        idx2 = i+ntmp;sidx2 = 2*ntmp

        idy1 = j-ntmp;sidy1 = 0
        idy2 = j+ntmp;sidy2 = 2*ntmp

        if ((idx1<0)|(idx2>=nnx)|(idy1<0)|(idy2>=nny)):
            break

        buf[sidx1:sidx2,sidy1:sidy2] = matrix[idx1:idx2,idy1:idy2]
        buf_x1[sidx1:sidx2,sidy1:sidy2] = xx1[idx1:idx2,idy1:idy2]
        buf_x2[sidx1:sidx2,sidy1:sidy2] = xx2[idx1:idx2,idy1:idy2]
        buf_tot = np.sum(buf)
        buf_bar = buf_tot/(4.0*ntmp*ntmp)

        buf_nrm = buf-buf_bar

        buf_nrm[buf_nrm <=0] = 0
        buf_nrm_tot = np.sum(buf_nrm)

        if buf_nrm_tot == 0.0 :
            break

        xc2 = np.sum(buf_nrm*buf_x1)/buf_nrm_tot
        xc1 = np.sum(buf_nrm*buf_x2)/buf_nrm_tot

        dx1 = x1-xc1
        dx2 = x2-xc2

        dr = np.sqrt(dx1*dx1+dx2*dx2)
        x1 = xc1
        x2 = xc2

        kk = kk+1

        if kk > 100 :
                break

    return x1,x2

#--------------------------------------------------------------------
def main():
    re = 1.0 # in units of arcsec
    boxsize = 6.0*re # in the units of Einstein Radius
    nnn = 1024
    dsx = boxsize/nnn

    xx01 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
    xx02 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
    xi2,xi1 = np.meshgrid(xx01,xx02)
    #----------------------------------------------------------------------
    g_amp = 1.0   	# peak brightness value
    g_sig = 0.02  	# Gaussian "sigma" (i.e., size)
    g_xcen = 0.03  	# x position of center (also try (0.0,0.14)
    g_ycen = 0.1  	# y position of center
    g_axrat = 1.0 	# minor-to-major axis ratio
    g_pa = 0.0    	# major-axis position angle (degrees) c.c.w. from x axis
    gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
    #----------------------------------------------------------------------
    #g_source = 0.0*xi1
    #g_source = gauss_2d(xi1,xi2,gpar) # modeling source as 2d Gaussian with input parameters.
    #----------------------------------------------------------------------
    xc1 = 0.0       #x coordinate of the center of lens (in units of Einstein radius).
    xc2 = 0.0       #y coordinate of the center of lens (in units of Einstein radius).
    q   = 0.7       #Ellipticity of lens.
    rc  = 0.1       #Core size of lens (in units of Einstein radius).
    re  = 1.0       #Einstein radius of lens.
    pha = 45.0      #Orintation of lens.
    lpar = np.asarray([xc1,xc2,q,rc,re,pha])
    #----------------------------------------------------------------------
    ai1,ai2,mua = lens_equation_sie(xi1,xi2,lpar)
    yi1 = xi1-ai1
    yi2 = xi2-ai2
    #----------------------------------------------------------------------


    gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
    g_lensimage = gauss_2d(yi1,yi2,gpar)
    #g_lensimage = gauss_2d(xi1,xi2,gpar)

    #----------------------------------------------------------------------
    g_amp = 5.0   	# peak brightness value
    g_sig = 0.5  	# Gaussian "sigma" (i.e., size)
    g_xcen = 0.03  	# x position of center (also try (0.0,0.14)
    g_ycen = 0.1  	# y position of center
    g_axrat = 1.0 	# minor-to-major axis ratio
    g_pa = 0.0    	# major-axis position angle (degrees) c.c.w. from x axis
    gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
    #----------------------------------------------------------------------
    g_source = gauss_2d(xi1,xi2,gpar) # modeling source as 2d Gaussian with input parameters.

    #add noise
    #g_noise = np.random.random_sample([nnn,nnn])
    g_noise = np.random.normal(0,1,[nnn,nnn])*1.0
    g_lensimage = g_lensimage+g_source+g_noise

    #smooth
    #g_lensimage = snf.uniform_filter(g_lensimage,size=4)
    #g_lensimage = snf.gaussian_filter(g_lensimage,1.0)

    levels = [0.0,1.0,2.0,3.0,4.0,5.0,6.0]
    pl.figure()
    pl.contourf(g_lensimage,levels)
    pl.colorbar()
    #----------------------------------------------------------------------
    nns = 16
    #dss = nns*dsx

    #x2 = 0.85
    #x1 = 0.095

    #xn1,xn2 = pixel_trans(x1,x2,xi1,xi2,g_lensimage,dss)

    xp01 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn/nns)+0.5*dsx*nns
    xp02 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn/nns)+0.5*dsx*nns
    xp2,xp1 = np.meshgrid(xp01,xp02)

    xp2 = xp2.reshape((nnn/nns*nnn/nns))
    xp1 = xp1.reshape((nnn/nns*nnn/nns))

    xr2 = xp2*0.0
    xr1 = xp1*0.0

    for i in xrange(len(xp2)):

        xr1[i],xr2[i] = pixel_trans(xp1[i],xp2[i],xi1,xi2,g_lensimage,nns)

    pl.figure()
    #pl.contourf(xi2,xi1,g_lensimage)
    pl.plot(xp1,xp2,'ro')

    pl.figure()
    #pl.contourf(xi2,xi1,g_lensimage)
    pl.plot(xr1,xr2,'ko')

    pl.show()

    #g_noise = np.random.random_sample([nnn,nnn])
    #g_noise = np.random.normal(0,1,[nnn,nnn])
    #print np.shape(g_noise)
    #g_lensimage = g_lensimage+g_noise

    #g_lensimage = snf.uniform_filter(g_lensimage,size=10)
    #g_lensimage = snf.gaussian_filter(g_lensimage,10.0)

    #gdif2,gdif1 = np.gradient(g_lensimage)

    #file_out = "./lensed_image.fits"
    #pyfits.writeto(file_out,g_lensimage,clobber=True)
    ##--------------------------lens images contour------------------------
    #levels = [0.15,0.30,0.45,0.60,0.75,0.9,1.05]
    #pl.figure(num=None,figsize=(15,5),dpi=80, facecolor='w', edgecolor='k')

    #a = pl.axes([0.05,0.1,0.26,0.8])
    #a.set_xlim(-boxsize/4.0,boxsize/4.0)
    #a.set_ylim(-boxsize/4.0,boxsize/4.0)
    #a.contourf(xi1,xi2,g_source,levels)
    #a.contour(yi1,yi2,mua,0,colors=('g'),linewidths = 2.0)

    #b = pl.axes([0.375,0.1,0.26,0.8])
    #b.set_xlim(-boxsize/4.0,boxsize/4.0)
    #b.set_ylim(-boxsize/4.0,boxsize/4.0)
    #b.contourf(xi1,xi2,g_lensimage,levels)
    #b.contour(xi1,xi2,mua,colors=('k'),linewidths = 2.0)
    ##savefig('output.eps')

    #c = pl.axes([0.70,0.1,0.26,0.8])
    #c.set_xlim(-boxsize/4.0,boxsize/4.0)
    #c.set_ylim(-boxsize/4.0,boxsize/4.0)
    #nn = 128
    #i = np.linspace(0,nn-1,nn)*nnn/nn
    #isub1,isub2 = np.meshgrid(i,i)
    #isub1 = np.array(isub1).astype(int)
    #isub2 = np.array(isub2).astype(int)

    #X = xi1[isub1,isub2]
    #Y = xi2[isub1,isub2]
    #U = gdif1[isub1,isub2]
    #V = gdif2[isub1,isub2]
    ##c.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=1)
    #c.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=0.3)
    #c.set_xlim([-1.5,1.5])
    #c.set_ylim([-1.5,1.5])

    #pl.show()
    return 0
#------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
