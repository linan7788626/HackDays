import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

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

def lensing_signals_sie(x1, x2, lpar):
    # lpar = np.asarray([l_xcen,l_ycen,l_re,l_rc,l_axrat,l_pa])
    # x coordinate of the center of lens (in units of Einstein radius).
    xc1 = lpar[0]
    # y coordinate of the center of lens (in units of Einstein radius).
    xc2 = lpar[1]
    re = lpar[2]   # Einstein radius of lens.
    rc = lpar[3]   # Core size of lens (in units of Einstein radius).
    q = lpar[4]   # Axis ratio of lens.
    pha = lpar[5]   # Orintation of lens.

    phirad = np.deg2rad(pha)
    cosa = np.cos(phirad)
    sina = np.sin(phirad)

    xt1 = (x1 - xc1) * cosa + (x2 - xc2) * sina
    xt2 = (x2 - xc2) * cosa - (x1 - xc1) * sina

    phi = np.sqrt(xt2 * xt2 + xt1 * q * xt1 * q + rc * rc)
    sq = np.sqrt(1.0 - q * q)
    pd1 = phi + rc / q
    pd2 = phi + rc * q
    fx1 = sq * xt1 / pd1
    fx2 = sq * xt2 / pd2
    qs = np.sqrt(q)

    a1 = qs / sq * np.arctan(fx1)
    a2 = qs / sq * np.arctanh(fx2)

    xt11 = cosa
    xt22 = cosa
    xt12 = sina
    xt21 = -sina

    fx11 = xt11 / pd1 - xt1 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd1 * pd1)
    fx22 = xt22 / pd2 - xt2 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd2 * pd2)
    fx12 = xt12 / pd1 - xt1 * \
        (xt1 * q * q * xt12 + xt2 * xt22) / (phi * pd1 * pd1)
    fx21 = xt21 / pd2 - xt2 * \
        (xt1 * q * q * xt11 + xt2 * xt21) / (phi * pd2 * pd2)

    a11 = qs / (1.0 + fx1 * fx1) * fx11
    a22 = qs / (1.0 - fx2 * fx2) * fx22
    a12 = qs / (1.0 + fx1 * fx1) * fx12
    a21 = qs / (1.0 - fx2 * fx2) * fx21

    rea11 = (a11 * cosa - a21 * sina) * re
    rea22 = (a22 * cosa + a12 * sina) * re
    rea12 = (a12 * cosa - a22 * sina) * re
    rea21 = (a21 * cosa + a11 * sina) * re

    #kappa = 0.5 * (rea11 + rea22)
    #shear1 = 0.5 * (rea12 + rea21)
    #shear2 = 0.5 * (rea11 - rea22)

    y11 = 1.0 - rea11
    y22 = 1.0 - rea22
    y12 = 0.0 - rea12
    y21 = 0.0 - rea21

    jacobian = y11 * y22 - y12 * y21
    mu = 1.0 / jacobian

    alpha1 = (a1 * cosa - a2 * sina) * re
    alpha2 = (a2 * cosa + a1 * sina) * re

    return alpha1, alpha2, mu#, kappa, shear1, shear2, mu

def lensed_images(x1, x2, lpar, gpar):
    al1, al2, mu = lensing_signals_sie(x1,x2,lpar)
    yi1 = x1-al1
    yi2 = x2-al2

    g_lensimage = 0.0*yi1
    g_lensimage = gauss_2d(yi1,yi2,gpar)

    return g_lensimage, mu, yi1, yi2
#--------------------------------------------------------------------
boxsize = 8.0 # (arcsec)
nnn = 128
dsx = boxsize/nnn

xi1 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
xi2 = np.linspace(-boxsize/2.0,boxsize/2.0,nnn)+0.5*dsx
xi1,xi2 = np.meshgrid(xi1,xi2)
#----------------------------------------------------------------------
l_xcen = 0.0  	# x position of center (also try (0.0,0.14)
l_ycen = 0.0  	# y position of center
l_re = 2.0   # Einstein radius of lens.
l_rc = 0.0   # Core size of lens (in units of Einstein radius).
l_axrat = 0.7   # Axis ratio of lens.
l_pa = 49.0   # Orintation of lens.

lpar = np.asarray([l_xcen,l_ycen,l_re,l_rc,l_axrat,l_pa])
#----------------------------------------------------------------------
g_amp = 1.0   	# peak brightness value
g_sig = 0.1  	# Gaussian "sigma" (i.e., size)
g_xcen = 0.0  	# x position of center (also try (0.0,0.14)
g_ycen = 0.0  	# y position of center
g_axrat = 0.7 	# minor-to-major axis ratio
g_pa = 45.0    	# major-axis position angle (degrees) c.c.w. from x axis

gpar = np.asarray([g_amp,g_sig,g_xcen,g_ycen,g_axrat,g_pa])
#----------------------------------------------------------------------
g_limage = 0.0*xi1
mua = 0.0*xi1
g_limage, mua,yi1,yi2 = lensed_images(xi1,xi2,lpar,gpar)

plt.figure(figsize=(6,6))
#fig, ax = plt.subplots()
#plt.subplots_adjust(left=0.25, bottom=0.25)
levels = [0.5,]
levels_mu = [0.0,]

plt.contour(xi1,xi2,g_limage,levels,colors = ('b'))
plt.contour(xi1,xi2,mua,levels_mu,colors = ('r'))
contour_axis = plt.gca()

plt.axis([-4, 4, -4, 4])

axcolor = 'lightgoldenrodyellow'
ax_re = plt.axes([0.25, 0.05, 0.65, 0.03], axisbg=axcolor)
ax_x1 = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
ax_x2 = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)

re_bar = Slider(ax_re, r"$R_e$",  0.0, 4.0, valinit=2.0)
x1_bar = Slider(ax_x1, r"$X_1$", -4.0, 4.0, valinit=0.0)
x2_bar = Slider(ax_x2, r"$X_2$", -4.0, 4.0, valinit=0.0)


def update(val):
    lpar[2] = re_bar.val
    gpar[2] = x1_bar.val
    gpar[3] = x2_bar.val

    contour_axis.clear()
    #contour_axis.contourf(xi1,xi2,lensed_images(xi1, xi2, lpar, gpar))
    g_tmp,mu_tmp,yi1_tmp,yi2_tmp = lensed_images(xi1, xi2, lpar, gpar)
    contour_axis.contour(xi1,xi2,g_tmp,levels,colors = ('b'))
    contour_axis.contour(xi1,xi2,mu_tmp,levels_mu,colors = ('r'))
    plt.draw()

re_bar.on_changed(update)
x1_bar.on_changed(update)
x2_bar.on_changed(update)

resetax = plt.axes([0.8, 0.0, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    x1_bar.reset()
    x2_bar.reset()
    re_bar.reset()
button.on_clicked(reset)

plt.show()
