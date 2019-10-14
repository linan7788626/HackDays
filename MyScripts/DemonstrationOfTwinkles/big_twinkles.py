#!/usr/bin/env python
import pygame
from pygame.locals import *
from sys import exit
import numpy as np

def nie_all(xi1,xi2,xc1,xc2,b,s,q,rot,ys1,ys2):

    x1,x2 = xy_rotate(xi1,xi2,xc1,xc2,rot)

    wx = np.sqrt(q*q*(x1*x1+s*s)+x2*x2)

    al1 = b/np.sqrt(1-q*q)*np.arctan(x1*np.sqrt(1-q*q)/(wx+s))
    al2 = b/np.sqrt(1-q*q)*np.arctanh(x2*np.sqrt(1-q*q)/(wx+q*q*s))

    kappa = b/(2.0*wx)

    hx = np.sqrt((wx+s)**2.0+(1-q*q)*x1*x1)
    phi = x1*al1+x2*al2-b*s*np.log(hx)+b*q*s*np.log((1+q)*s)

    Kc = 1.0
    #Kc = (1.0+zl)/c*(Dl*Ds/Dls)
    td = Kc*(0.5*((al1)**2.0+(al2)**2.0)-phi)
    #td = Kc*(0.5*((x1-ys1)**2.0+(x2-ys2)**2.0)-phi)

    y1 = x1-al1
    y2 = x2-al2

    y1,y2 = xy_rotate(y1,y2,xc1,xc2,-rot)

    demon1 = ((wx+s)**2+(1.0-q*q)*x1*x1)*wx
    demon2 = (((wx+q*q*s)**2-(1.0-q*q)*x2*x2)*wx)
    y11 = 1-b*(wx*(wx+s)-q*q*x1*x1)/demon1
    y22 = 1-b*(wx*(wx+q*q*s)-x2*x2)/demon2
    y12 = -b*x1*x2/demon1
    y21 = -b*x1*x2*q*q/demon2

    mu = 1.0/(y11*y22-y12*y21)

    return phi,td,al1,al2,kappa,mu,y1,y2

def lensed_images(xi1,xi2,yi1,yi2,gpar):
    g_image = gauss_2d(xi1,xi2,gpar)
    g_lensimage = gauss_2d(yi1,yi2,gpar)
    return g_image,g_lensimage

def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
    return (xnew,ynew)

def gauss_2d(x, y, par):
    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / np.abs(par[1])**2
    return par[0] * np.exp(-0.5*r_ell_sq)

def parabola_1d(x,xb,xc,a):
    res = a*(x-xb)*(2.0*xc-(x-xb))
    res[res<=0] = 0.0
    return res

def lens_galaxies(xi1,xi2,glpar):
    g_lens = gauss_2d(xi1,xi2,glpar)
    return g_lens

def main():

    nnn = 512
    boxsize = 4.0
    dsx = boxsize/nnn
    xi1 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi2 = np.linspace(-boxsize/2.0,boxsize/2.0-dsx,nnn)+0.5*dsx
    xi1,xi2 = np.meshgrid(xi1,xi2)

    pygame.init()
    FPS = 60
    fpsClock = pygame.time.Clock()
    screen = pygame.display.set_mode((nnn, nnn), 0, 32)
    pygame.display.set_caption("Gravitational Lensing Toy")
    mouse_cursor = pygame.Surface((nnn,nnn))

    #----------------------------------------------------

    base0 = np.zeros((nnn,nnn,3),'uint8')
    base1 = np.zeros((nnn,nnn,3),'uint8')
    base2 = np.zeros((nnn,nnn,3),'uint8')

    #----------------------------------------------------
    # lens parameters for main halo
    xlc1 = 0.0
    xlc2 = 0.0
    ql0 = 0.699999999999
    rc0 = 0.000000000001
    re0 = 1.0
    phi0 = 30.0
    lpar = np.asarray([xlc1, xlc2, re0, rc0, ql0, phi0])

    lpars_list = []
    lpars_list.append(lpar)
    #----------------------------------------------------
    # lens parameters for main halo
    xls1 = 0.7
    xls2 = 0.8
    qls = 0.999999999999
    rcs = 0.000000000001
    res = 0.0
    phis = 0.0
    lpars = np.asarray([xls1, xls2, res, rcs, qls, phis])
    lpars_list.append(lpars)

    ap0 = 1.0
    l_sig0 = 0.5
    glpar  = np.asarray([ap0,l_sig0,xlc1,xlc2,ql0,phi0])

    g_lens = lens_galaxies(xi1,xi2,glpar)

    base0[:,:,0] = g_lens*256
    base0[:,:,1] = g_lens*128
    base0[:,:,2] = g_lens*0

    x = 0.33984375*nnn/2
    y = -0.11328125*nnn/2
    step = 1
    gr_sig = 0.1

    LeftButton=0

    #----------------------------------------------------

    ic = FPS/6.0

    i = 0
    while True:
        i = i+1
        for event in pygame.event.get():
            if event.type == QUIT:
                exit()
            if event.type == MOUSEMOTION:

                if event.buttons[LeftButton]:
                    rel = event.rel
                    x += rel[0]
                    y += rel[1]

            #----------------------------------------------
            if event.type == pygame.MOUSEBUTTONDOWN:
                if event.button == 4:
                    gr_sig -= 0.1
                    if gr_sig <0.01:
                        gr_sig = 0.01

                elif event.button == 5:
                    gr_sig += 0.01
                    if gr_sig >0.4:
                        gr_sig = 0.4



        keys = pygame.key.get_pressed()  #checking pressed keys
        if keys[pygame.K_RIGHT]:
            x += step
            if x > 500:
                x = 500
        if keys[pygame.K_LSHIFT] & keys[pygame.K_RIGHT]:
            x += 30*step

        if keys[pygame.K_LEFT]:
            x -= step
            if x < -500:
                x = -500

        if keys[pygame.K_LSHIFT] & keys[pygame.K_LEFT]:
            x -= 30*step

        if keys[pygame.K_UP]:
            y -= step
            if y < -500 :
                y = -500
        if keys[pygame.K_LSHIFT] & keys[pygame.K_UP]:
            y -= 30*step

        if keys[pygame.K_DOWN]:
            y += step
            if y > 500 :
                y = 500
        if keys[pygame.K_LSHIFT] & keys[pygame.K_DOWN]:
            y += 30*step


        #----------------------------------------------
        if keys[pygame.K_MINUS]:
            gr_sig -= 0.01
            if gr_sig <0.01:
                gr_sig = 0.01

        if keys[pygame.K_EQUALS]:
            gr_sig += 0.01
            if gr_sig >0.1:
                gr_sig = 0.1

        #gr_sig = 0.005

        #----------------------------------------------
        #parameters of source galaxies.
        #----------------------------------------------
        g_amp = 1.0         # peak brightness value
        g_sig = gr_sig          # Gaussian "sigma" (i.e., size)
        g_xcen = x*2.0/nnn  # x position of center
        g_ycen = y*2.0/nnn  # y position of center
        g_axrat = 1.0       # minor-to-major axis ratio
        g_pa = 0.0          # major-axis position angle (degrees) c.c.w. from y axis
        gpar = np.asarray([g_amp, g_sig, g_ycen, g_xcen, g_axrat, g_pa])
        #----------------------------------------------

        #----------------------------------------------
        #parameters of SNs.
        #----------------------------------------------
        g_amp = 1.0         # peak brightness value
        g_sig = 0.05          # Gaussian "sigma" (i.e., size)
        g_xcen = x*2.0/nnn+0.05  # x position of center
        g_ycen = y*2.0/nnn+0.05  # y position of center
        g_axrat = 1.0       # minor-to-major axis ratio
        g_pa = 0.0          # major-axis position angle (degrees) c.c.w. from y axis
        gpsn = np.asarray([g_amp, g_sig, g_ycen, g_xcen, g_axrat, g_pa])

        phi,td,ai1,ai2,kappa,mu,yi1,yi2 = nie_all(xi1,xi2,xlc1,xlc2,re0,rc0,ql0,phi0,g_ycen,g_xcen)
        g_image,g_lensimage = lensed_images(xi1,xi2,yi1,yi2,gpar)
        g_image = g_image*0.0
        g_lensimage = g_lensimage*0.0
        g_sn,g_lsn = lensed_images(xi1,xi2,yi1,yi2,gpsn)

        sktd = td/td.max()*ic
        itmp = (i)%(FPS)
        ratio = parabola_1d(itmp,40+sktd,ic,2.0/ic**2.0)
        ratio0 = parabola_1d(itmp,0.0*sktd,ic,2.0/ic**2.0)

        base1[:,:,0] = g_sn*168*(1.0+ratio0)/2+g_image*256
        base1[:,:,1] = g_sn*168*(1.0+ratio0)/2+g_image*256
        base1[:,:,2] = g_sn*168*(1.0+ratio0)/2+g_image*256

        base2[:,:,0] = g_lsn*168*(1.0+ratio)/2+g_lensimage*102
        base2[:,:,1] = g_lsn*168*(1.0+ratio)/2+g_lensimage*178
        base2[:,:,2] = g_lsn*168*(1.0+ratio)/2+g_lensimage*256

        wf = base1+base2

        idx1 = wf>=base0
        idx2 = wf<base0

        base = base0*0
        base[idx1] = wf[idx1]
        base[idx2] = base0[idx2]

        pygame.surfarray.blit_array(mouse_cursor,base)
        screen.blit(mouse_cursor, (0, 0))
        pygame.display.update()
        fpsClock.tick(FPS)

if __name__ == '__main__':
    main()

