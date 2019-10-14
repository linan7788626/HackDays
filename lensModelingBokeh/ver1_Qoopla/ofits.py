import scipy.optimize as opt
import numpy as np

def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()

def lensing_signals_sie(x1, x2, xc1, xc2, re, rc, q, phl):

    phirad = np.deg2rad(phl)
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

    alpha1 = (a1 * cosa - a2 * sina) * re
    alpha2 = (a2 * cosa + a1 * sina) * re

    return alpha1, alpha2 #, kappa, shear1, shear2, mu

def create_images((x1,x2),xc1,xc2,re,rc,ql,phl,yc1,yc2,sigs,amp,qs,phs):

    if qs > 1.0 :
        qs = 1.0/qs
        phs = phs + 90

    if ql > 1.0 :
        ql = 1.0/ql
        phl = phl + 90

    sig1 = sigs#/qs#*0.693
    sig2 = sigs*qs#*0.693
    al1, al2 = lensing_signals_sie(x1, x2, xc1, xc2, re, 0.0, ql, phl)

    phs = np.deg2rad(phs)

    y1 = x1 - al1
    y2 = x2 - al2
    res = twoD_Gaussian((y1, y2), 1.0, yc1, yc2, sig1, sig2, phs, 0.0)
    return res

def optimize_pars(x1, x2, data, lpar_in, spar_in):
    # Create x and y indices
    #nnn = int(np.sqrt(len(data)))

    initial_guess = (lpar_in[0],lpar_in[1],lpar_in[2],lpar_in[3],lpar_in[4],lpar_in[5],spar_in[0],spar_in[1],spar_in[2],spar_in[3],spar_in[4],spar_in[5])

    popt, pcov = opt.curve_fit(create_images, (x1, x2), data, p0=initial_guess)
    lpar_new = popt[:6]
    spar_new = popt[6:]

    return lpar_new, spar_new

#if __name__ == '__main__':
    ## Create x and y indices
    #nnn = 400
    #dsx = 0.0225
    #bsz = dsx*nnn
    #x1 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
    #x2 = np.linspace(-bsz/2.0, bsz/2.0-dsx, nnn)+dsx/2.0
    #x1, x2 = np.meshgrid(x1, x2)

    #fname = "./ering.jpg"
    #data = loadin_images(fname)


    #fname_json = sys.argv[1]
    #initial_guess_0 = load_guessed_pars(fname_json)

    #lpar_in = [initial_guess_0[0],initial_guess_0[1],initial_guess_0[2],0.00,initial_guess_0[3],initial_guess_0[4]]
    #spar_in = [initial_guess_0[5],initial_guess_0[6],initial_guess_0[7],1.00,initial_guess_0[8],initial_guess_0[9]]

    #initial_guess = [initial_guess_0[0],initial_guess_0[1],initial_guess_0[2],0.00,initial_guess_0[3],initial_guess_0[4],initial_guess_0[5],initial_guess_0[6],initial_guess_0[7],1.00,initial_guess_0[8],initial_guess_0[9]]


    #lpar_new, spar_new = optimize_pars(x1, x2, data, lpar_in, spar_in)

    #popt, pcov = opt.curve_fit(create_images, (x1, x2), data, p0=initial_guess)

    #lpar_new = popt[:6]
    #spar_new = popt[6:]

    #print type(lpar_new), lpar_new
    #print type(spar_new), spar_new

    #data_guessed = create_images((x1, x2), *initial_guess)

    #levels = [0.5]
    #fig, ax = plt.subplots(1, 1)
    #ax.hold(True)
    #ax.imshow(data.reshape(nnn, nnn), cmap=plt.cm.winter, origin='bottom',
        #extent=(x1.min(), x1.max(), x2.min(), x2.max()))
    #ax.contour(x1, x2, data_guessed.reshape(nnn, nnn), levels, colors='r')

    #data_fitted = create_images((x1, x2), *popt)

    #save_improved_pars(fname_json, popt, "Improved_"+fname_json)

    #fig, ax = plt.subplots(1, 1)
    #ax.hold(True)
    #ax.imshow(data.reshape(nnn, nnn), cmap=plt.cm.winter, origin='bottom',
        #extent=(x1.min(), x1.max(), x2.min(), x2.max()))
    #ax.contour(x1, x2, data_fitted.reshape(nnn, nnn), levels, colors='r')

    #plt.figure()
    #plt.contourf(x1,x2,(data_fitted-data).reshape(nnn, nnn), cmap=plt.cm.jet)
    #plt.colorbar()
    #plt.show()
