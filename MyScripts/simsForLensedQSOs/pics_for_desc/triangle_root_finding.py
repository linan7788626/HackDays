#!/usr/bin/env python
import numpy as np
import pylab as pl

#from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
#import scipy.optimize as sco

def build_triangles(xgrids1,xgrids2):

    indices_matrix = np.indices(np.shape(xgrids1))

    idx_vertex1 = indices_matrix[:,:-1,:-1]
    idx_vertex2 = indices_matrix[:,1:,:-1]
    idx_vertex3 = indices_matrix[:,1:,1:]
    idx_vertex4 = indices_matrix[:,:-1,1:]

    tri1 = [idx_vertex1,idx_vertex2,idx_vertex3]
    tri2 = [idx_vertex1,idx_vertex3,idx_vertex4]

    return np.array(tri1),np.array(tri2)

def sign_tri(p1,p2,p3):
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])

def PointInTriangle(pt,v1,v2,v3):
    b1 = sign_tri(pt, v1, v2) <= 0.0
    b2 = sign_tri(pt, v2, v3) <= 0.0
    b3 = sign_tri(pt, v3, v1) <= 0.0

    return ((b1 == b2) & (b2 == b3))

def triangle_vectors(ind,jnd,tri1,xgrids1,xgrids2):
    v1 = np.array([xgrids1[tri1[0,0,ind,jnd],tri1[0,1,ind,jnd]],xgrids2[tri1[0,0,ind,jnd],tri1[0,1,ind,jnd]]])
    v2 = np.array([xgrids1[tri1[1,0,ind,jnd],tri1[1,1,ind,jnd]],xgrids2[tri1[1,0,ind,jnd],tri1[1,1,ind,jnd]]])
    v3 = np.array([xgrids1[tri1[2,0,ind,jnd],tri1[2,1,ind,jnd]],xgrids2[tri1[2,0,ind,jnd],tri1[2,1,ind,jnd]]])

    return v1,v2,v3

def tria_vectors(i,j,xgrids1,xgrids2):

    i
    ip1 = i+1
    jp1 = j+1

    v1 = np.array([xgrids1[i,j],xgrids2[i,j]])
    v2 = np.array([xgrids1[ip1,j],xgrids2[ip1,j]])
    v3 = np.array([xgrids1[ip1,jp1],xgrids2[ip1,jp1]])
    return v1,v2,v3

def trib_vectors(i,j,xgrids1,xgrids2):

    i
    ip1 = i+1
    jp1 = j+1

    v1 = np.array([xgrids1[i,j],xgrids2[i,j]])
    v2 = np.array([xgrids1[ip1,jp1],xgrids2[ip1,jp1]])
    v3 = np.array([xgrids1[i,jp1],xgrids2[i,jp1]])
    return v1,v2,v3

def standard_grids(x0,y0,bsz,ngrids):

    dsx = bsz/ngrids
    xi1 = np.linspace(-bsz/2.0,bsz/2.0-dsx,ngrids)+0.5*dsx+x0
    xi2 = np.linspace(-bsz/2.0,bsz/2.0-dsx,ngrids)+0.5*dsx+y0
    xi1,xi2 = np.meshgrid(xi1,xi2)

    return xi1,xi2

def new_grids(x0,y0,bsz,ngrids):

    #dsx = bsz/ngrids
    #xi1 = np.linspace(-bsz/2.0,bsz/2.0-dsx,ngrids)+0.2*dsx+x0
    #xi2 = np.linspace(-bsz/2.0,bsz/2.0-dsx,ngrids)+0.2*dsx+y0
    xi1 = np.linspace(-bsz/2.0,bsz+bsz/2.0,ngrids)+x0
    xi2 = np.linspace(-bsz/2.0,bsz+bsz/2.0,ngrids)+y0
    #xi1 = np.linspace(x0,xf0,ngrids)
    #xi2 = np.linspace(y0,yf0,ngrids)
    xi1,xi2 = np.meshgrid(xi1,xi2)

    return xi1,xi2,xi1[1,1]-xi1[0,0]


def mapping_triangles(ys1,ys2,xgrids1,xgrids2,lgrids1,lgrids2):

    tri1,tri2 = build_triangles(xgrids1,xgrids2)
    ntris = np.shape(xgrids1)[0]-1
    #dsx = xgrids1[1,1]-xgrids1[0,0]
    xroot1 = []
    xroot2 = []

    ncount = 0
    for i in xrange(ntris):
        for j in xrange(ntris):
            xv1,xv2,xv3 = tria_vectors(i,j,xgrids1,xgrids2)
            lv1,lv2,lv3 = tria_vectors(i,j,lgrids1,lgrids2)
            if PointInTriangle([ys1,ys2],lv1,lv2,lv3):
                xroot1.append(xv1[0])
                xroot2.append(xv1[1])
                ncount = ncount + 1

                #pl.figure(figsize=(10,10))
                #pl.xlim(-2.0,2.0)
                #pl.ylim(-2.0,2.0)
                #pl.plot([xv1[0],xv2[0],xv3[0]],[xv1[1],xv2[1],xv3[1]],'ro')
                #pl.plot([lv1[0],lv2[0],lv3[0]],[lv1[1],lv2[1],lv3[1]],'go')
                #pl.plot(ys1,ys2,'ko')

    for i in xrange(ntris):
        for j in xrange(ntris):
            xv1,xv2,xv3 = trib_vectors(i,j,xgrids1,xgrids2)
            lv1,lv2,lv3 = trib_vectors(i,j,lgrids1,lgrids2)
            if PointInTriangle([ys1,ys2],lv1,lv2,lv3):
                #xroot1.append((xv1[0]+xv2[0]+xv3[0])/3.0)
                #xroot2.append((xv1[1]+xv2[1]+xv3[1])/3.0)
                xroot1.append(xv1[0])#+dsx/2.0)
                xroot2.append(xv1[1])#+dsx/2.0)
                ncount = ncount + 1

                #print xv2[0]-xv1[0],xv3[0]-xv1[0],dsx
                #print xv2[1]-xv1[1],xv3[1]-xv1[1],dsx

                #pl.figure(figsize=(10,10))
                #pl.xlim(-2.0,2.0)
                #pl.ylim(-2.0,2.0)
                #pl.plot([xv1[0],xv2[0],xv3[0]],[xv1[1],xv2[1],xv3[1]],'ro')
                #pl.plot([lv1[0],lv2[0],lv3[0]],[lv1[1],lv2[1],lv3[1]],'go')
                #pl.plot(ys1,ys2,'ko')



    return np.array(xroot1),np.array(xroot2),ncount

def isNewImage(x1,x2,xil1,xil2):
    rdist = np.hypot((x1-xil1),(x2-xil2))
    lidx = len(rdist[rdist < 1e-7])
    return lidx

def detect_local_maxima(image):
    neighborhood = generate_binary_structure(2,2)
    local_max = maximum_filter(image, footprint=neighborhood)==image
    background = (image==0)
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
    detected_peaks = local_max - eroded_background
    return detected_peaks

def hfunc(x1,x2,rcore,qe):
    res = np.sqrt(qe*qe*(rcore*rcore+x1*x1)+x2*x2)
    return res


def nie_phi(x1,x2,re0,rcore,qe):

    res0 = re0/np.sqrt(1-qe*qe)
    al1 = res0*np.arctan(x1*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore))
    al2 = res0*np.arctanh(x2*np.sqrt(1-qe*qe)/(hfunc(x1,x2)+rcore*qe*qe))

    res1 = x1*al1+x2*al2
    res2 = re0*rcore*np.log(np.sqrt((hfunc(x1,x2)+rcore)**2.0+(1-qe*qe)*x1*x1))
    res = res1-res2
    return res

def nie_alphas(xx1,xx2,xc1,xc2,re0,rcore,qe):
    x1 = xx1-xc1
    x2 = xx2-xc2
    res0 = re0/np.sqrt(1-qe*qe)
    al1 = np.arctan(x1*np.sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore))
    al2 = np.arctanh(x2*np.sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore*qe*qe))
    return res0*al1,res0*al2

def nie_kappa(x1,x2,re0,rcore,qe):
    res = re0/(2.0*np.sqrt(qe*qe*(rcore+x1*x1)+x2*x2))
    return res

def nie_mu(x1,x2,re0,rcore,qe):
    res = 1.0/(1.0-re0/hfunc(x1,x2)-re0*re0*rcore/(hfunc(x1,x2)*((hfunc(x1,x2)+rcore)**2+(1-qe*qe)*x1*x1)))
    return res

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
    td2 = Kc*(0.5*((x1-ys1)**2.0+(x2-ys2)**2.0)-phi)

    y1 = x1-al1
    y2 = x2-al2

    y1,y2 = xy_rotate(y1,y2,xc1,xc2,-rot)

#------------------------------------------------------------------
    demon1 = ((wx+s)**2+(1.0-q*q)*x1*x1)*wx
    demon2 = (((wx+q*q*s)**2-(1.0-q*q)*x2*x2)*wx)
    y11 = 1-b*(wx*(wx+s)-q*q*x1*x1)/demon1
    y22 = 1-b*(wx*(wx+q*q*s)-x2*x2)/demon2
    y12 = -b*x1*x2/demon1
    y21 = -b*x1*x2*q*q/demon2

    mu = 1.0/(y11*y22-y12*y21)

    yss1 = y1-ys1
    yss2 = y2-ys2

    return phi,td,al1,al2,kappa,mu,y1,y2,td2,yss1,yss2

def lensed_images(xi1,xi2,yi1,yi2,gpar):

    g_image = gauss_2d(xi1,xi2,gpar)
    g_lensimage = gauss_2d(yi1,yi2,gpar)

    return g_image,g_lensimage

def lensed_images_point(xi1,xi2,yi1,yi2,gpar):
    g_image = tophat_2d(xi1,xi2,gpar)
    g_lensimage = tophat_2d(yi1,yi2,gpar)
    return g_image,g_lensimage

def xy_rotate(x, y, xcen, ycen, phi):

    phirad = np.deg2rad(phi)
    xnew = (x - xcen) * np.cos(phirad) + (y - ycen) * np.sin(phirad)
    ynew = (y - ycen) * np.cos(phirad) - (x - xcen) * np.sin(phirad)
    return (xnew,ynew)

def tophat_2d(x, y, par):

    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4])/np.abs(par[1])**2
    res = r_ell_sq*0.0
    res[r_ell_sq<=1.0] = par[0]
    return res

def gauss_2d(x, y, par):

    (xnew,ynew) = xy_rotate(x, y, par[2], par[3], par[5])
    r_ell_sq = ((xnew**2)*par[4] + (ynew**2)/par[4]) / np.abs(par[1])**2
    return par[0] * np.exp(-0.5*r_ell_sq)

def gauss_1d(x, x0,sigma,a):

    r_ell_sq = ((x-x0)**2.0/sigma**2)
    res = a*np.exp(-0.5*r_ell_sq)
    return res

def parabola_1d(x,xb,xc,a):
    res = a*(x-xb)*(2.0*xc-(x-xb))
    res[res<=0] = 0.0
    return res

def find_critical_curve(mu):
    rows,cols = np.indices(np.shape(mu))
    cdtn = np.sign(mu)*(np.sign(mu[rows-1,cols])+np.sign(mu[rows,cols-1])+np.sign(mu[(rows+1)%len(rows),cols])+np.sign(mu[rows,(cols+1)%len(cols)]))

    res = mu*0
    res[cdtn<4] = 1
    res[cdtn>=4] = 0

    return res

def lens_galaxies(xi1,xi2,glpar):

    g_lens = gauss_2d(xi1,xi2,glpar)

    return g_lens

def roots_zeros(xi1,xi2,ai1,ai2,ys1,ys2):

    yi1 = xi1-ai1
    yi2 = xi2-ai2

    xroot1,xroot2,nroots = mapping_triangles(ys1,ys2,xi1,xi2,yi1,yi2)

    return xroot1,xroot2,nroots

def refine_roots(xb1,xb2,xl1,xl2,bsz,nnn,ys1,ys2):

    ql0 = 0.699999999999
    rc0 = 0.000000000000
    re0 = 1.0

    xi1,xi2,dsx = new_grids(xb1,xb2,bsz,nnn)
    alpha1,alpha2 = nie_alphas(xi1,xi2,xl1,xl2,re0,rc0,ql0)
    yi1 = xi1-alpha1
    yi2 = xi2-alpha2


    xroot1,xroot2,nroots = mapping_triangles(ys1,ys2,xi1,xi2,yi1,yi2)

    #pl.figure(figsize=(10,10))
    #pl.xlim(-2.0,2.0)
    #pl.ylim(-2.0,2.0)
    #pl.plot(yi1,yi2,'ro')
    #pl.plot(xi1,xi2,'go')
    #pl.plot(ys1,ys2,'ko')

    return xroot1,xroot2,nroots,dsx
def run_main():
    nnn = 128
    bsz = 4.0

    xb1 = 0.0
    xb2 = 0.0

    xl1 = 0.0
    xl2 = 0.0

    re0 = 1.0
    rc0 = 0.0
    ql0 = 0.7

    ys1 = np.random.random(1)*(0.2-0.1)+0.05
    ys2 = np.random.random(1)*(0.2-0.1)+0.05

    xi1,xi2 = standard_grids(xb1,xb2,bsz,nnn)
    ai1,ai2 = nie_alphas(xi1,xi2,xl1,xl2,re0,rc0,ql0)

    xroot1,xroot2,nroots = roots_zeros(xi1,xi2,ai1,ai2,ys1,ys2)

    print xroot1,xroot2
    pl.figure(figsize=(10,10))
    pl.xlim(-2.0,2.0)
    pl.ylim(-2.0,2.0)
    pl.plot(xroot1,xroot2,'go')
    pl.plot(ys1,ys2,'ko')
    return 0

if __name__ == '__main__':
    run_main()
    pl.show()
