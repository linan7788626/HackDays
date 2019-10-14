import numpy as np
import scipy.interpolate
import scipy.ndimage
import scipy.signal
import scipy.misc
import scipy.ndimage.filters
import pyfits
import congrid
from pylab import *

#-----------------------------------------------------------------------
def plotim(imin):
    im=imin.flat
    idx=np.random.uniform(0,1,100000)*len(im)
    idx=idx.astype(int)
    idx1=np.abs(im[idx]) > 1e-4

    idx=idx[idx1]
    a = im[idx]
    a.sort()
    im = a
    im=im[0.2*len(im):0.8*len(im)]

    sky=np.median(im)
    skysig=np.std(im)

    #print sky, skysig

    #sky = np.median(im)
    #skysig = np.std(im)

    #bin=2
    #bin1=sky-skysig*10
    #bin2=sky+skysig*10
    #bin=skysig/10.0
    #h=histogram(im[idx],min=bin1,max=bin2,binsize=bin)
    #hx=findgen(len(h))*bin+bin1+bin/2.0
    #param=[max(h),sky,skysig/0.82]
    #gf=gaussfit(hx,h,fit,estimates=param,nterms=3)

    #!p.multi=[0,1,2]
    #plot,hx,h,psym=10
    #oplot,hx,gf,col=255
    #foo=where(hx lt fit(1)+fit(2))
    #gf1=gaussfit(hx(foo),h(foo),fit,nterms=3,estimates=param)
    #oplot,hx,gf1,col='00ff00'x
    #foo1=where(hx lt fit(1)+fit(2)*0.25)
    #gf2=gaussfit(hx(foo1),h(foo1),fit,nterms=3,estimates=param)
    #oplot,hx,gf2,col='00ffff'x

    #plot,hx,(h-gf)/max(h),col=255
    #oplot,hx(foo),(h(foo)-gf1)/max(h),col='00ff00'x
    #oplot,hx(foo1),(h(foo1)-gf2)/max(h),col='00ffff'x

    #oplot,[0,0]+fit(1),[-100,100]
    #sky=fit(1)
    #skysig=fit(2)
    return sky,skysig

#-----------------------------------------------------------------------
def grab_individual_bright_stars():

    sptb=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.g.fits')
    sptv=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.r.fits')
    spti=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.i.fits')

    #0307-5042 center :  5441 5372 in r: g,i otherwise
    #observed plate scale is 0.16 "/pixel
    #PISCO planned scale is 0.157 "/pixel
    #plate scale very similar to REBINNED TO PISCO simulation, so need
    #1174x1174 cutout from the real MEgaCam data

    #STAR1, box +-110
    #1814,6063 in I
    #2190,6173 in V
    #1815,6063 in B
    star1i=spti[6063-110:6063+110,1814-110:1814+110]
    star1v=sptv[6173-110:6173+110,2190-110:2190+110]
    star1b=sptb[6063-110:6063+110,1815-110:1815+110]

    war1=np.zeros((221+30,221+30))
    war1[15:15+220,15:15+220]=1
    war1=scipy.ndimage.filters.uniform_filter(war1,15)
    war1=war1[15:15+220,15:15+220]
    nar1=1-war1

    #STAR2, box +-45
    #6354,7092 in I
    #6730,7202 in V
    #6355,7092 in B
    star2i=spti[7092-45:7092+45,6354-45:6354+45]
    star2v=sptv[7202-45:7202+45,6730-45:6730+45]
    star2b=sptb[7092-45:7092+45,6355-45:6355+45]

    war2=np.zeros((91+30,91+30))
    war2[15:15+90,15:15+90]=1
    war2=scipy.ndimage.filters.uniform_filter(war2,15)
    war2=war2[15:15+90,15:15+90]
    nar2=1-war2

    #STAR3, box +-43
    #2531,3440 in I
    #2907,3550 in V
    #2532,3440 in B

    star3i=spti[3440-43:3440+43,2531-43:2531+43]
    star3v=sptv[3550-43:3550+43,2907-43:2907+43]
    star3b=sptb[3440-43:3440+43,2532-43:2532+43]

    war3=np.zeros((87+30,87+30))
    war3[15:15+86,15:15+86]=1
    war3=scipy.ndimage.filters.uniform_filter(war3,15)
    war3=war3[15:15+86,15:15+86]
    nar3=1-war3

    return star1i,star1v,star1b,war1,nar1,star2i,star2v,star2b,war2,nar2,star3i,star3v,star3b,war3,nar3



    #---------grab individual bright stars from real data--------

    #sptb=sptb(5066-587:5066+586,5262-587:5262+586)
    #sptv=sptv(5441-587:5441+586,5372-587:5372+586)
    #spti=spti(5065-587:5065+586,5262-587:5262+586)

    #plotim,sptb,skypb,sigpb
    #plotim,sptv,skypv,sigpv
    #plotim,spti,skypi,sigpi

    #r-band moffat fit 3.0 pix
    #g-band moffat fit 5.5 pix
    #i-band moffat fit 4.5 pix

    #in rebinned pix, the moffat
    #0.5 psf = 3.2 pix
    #0.6 psf = 3.8 pix
    #0.7 psf = 4.5 pix
    #0.8 psf = 5.1 pix
    #0.9 psf = 5.7 pix

def read_PSF_kernels():

    #psfs are done as for 0.09" pixels
    psfv=pyfits.getdata('psf_moff_0.5.fits')
    psfv=psfv/np.sum(psfv)
    psfi=pyfits.getdata('psf_moff_0.7.fits')
    psfi=psfi/np.sum(psfi)
    psf1=pyfits.getdata('psf_moff_0.8.fits')
    psf2=pyfits.getdata('psf_moff_0.9.fits')
    psfb=psf1+psf2
    psfb=psfb/np.sum(psfb)

    return psf1,psf2,psfi,psfv,psfb


#code top take Nan's input lensed images and make them as
#though they were observed by gemini...

#fwhm is 2.355x gaussian sigma
#PISCO 2x2 binned pixel 0.157"
#0.7" FWHM=4.46 pix
def read_mem_gals():
    #read cluster galaxy images from Lindsey's code
    tag = '1'
    clb=pyfits.getdata('testgalsim'+tag+'_g.fits.txt')
    clv=pyfits.getdata('testgalsim'+tag+'_r.fits.txt')
    cli=pyfits.getdata('testgalsim'+tag+'_i.fits.txt')
    clb=clb[3333-1024:3333+1024,3333-1024:3333+1024]
    clv=clv[3333-1024:3333+1024,3333-1024:3333+1024]
    cli=cli[3333-1024:3333+1024,3333-1024:3333+1024]

    return cli,clv,clb

def read_lensed_images():

    iim=pyfits.getdata('./bdded_540_I.fits')
    vim=pyfits.getdata('./bdded_540_V.fits')
    bim=pyfits.getdata('./bdded_540_B.fits')

    #scale lensed image pixels
    bim=bim*500.*1.3
    vim=vim*700.*1.3
    iim=iim*600.*1.3

    return iim,vim,bim

def read_los_images():

    bom=pyfits.getdata('./final_images/los_0_B.fits')
    vom=pyfits.getdata('./final_images/los_0_V.fits')
    iom=pyfits.getdata('./final_images/los_0_I.fits')

    #scale lensed image pixels
    bom=bom*500.*1.3
    vom=vom*700.*1.3
    iom=iom*600.*1.3

    return iom,vom,bom

def step_through_images(iim,vim,bim,iom,vom,bom,cli,clv,clb,star1i,star1v,star1b,star2i,star2v,star2b,star3i,star3v,star3b,war1,war2,war3,psfi,psfv,psfb):

    #scaling variables to add 'real' bright stars...
    war1=war1*0.25
    nar1=1-war1
    war2=war2*0.25
    nar2=1-war2
    war3=war3*0.25
    nar3=1-war3

    #in sim cluster gals
    #vflux = 0.440 iflux
    #bflux = 0.105 iflux
    #in SPT obs of z=0.55 cluster in gri
    #vflux = 1.09  iflux
    #bflux = 0.054 iflux
    #print np.shape(bom),np.shape(bim),np.shape(clb)

    #hence adjust summed gri frames to match spt fluxes << ***
    #here are adding simmed cluster gals with scaling, to lensed images
    skyi = 2.5515003687
    sigi = 6.69451075049
    skyv = 3.43434568331
    sigv = 8.14642639767
    skyb = 0.661582802206
    sigb = 2.38341579459

    save_to_jpeg(iim*5,vim*5,bim*5,2048,skyi,sigi,skyv,sigv,skyb,sigb,"slides_0.jpg")

    igm=iim+cli*1.00
    vgm=vim+clv*2.48
    bgm=bim+clb*0.52

    save_to_jpeg(igm*5,vgm*5,bgm*5,2048,skyi,sigi,skyv,sigv,skyb,sigb,"slides_1.jpg")

    #igm=igm+iom[2048-1024:2048+1024,2048-1024:2048+1024]
    #vgm=vgm+vom[2048-1024:2048+1024,2048-1024:2048+1024]
    #bgm=bgm+bom[2048-1024:2048+1024,2048-1024:2048+1024]

    #save_to_jpeg(igm*5,vgm*5,bgm*5,2048,skyi,sigi,skyv,sigv,skyb,sigb,"slides_2.jpg")

    #bgm=bim[2048-1024:2048+1024,2048-1024:2048+1024]+clb*0.52
    #vgm=vim[2048-1024:2048+1024,2048-1024:2048+1024]+clv*2.48
    #igm=iim[2048-1024:2048+1024,2048-1024:2048+1024]+cli*1.00

    #bgm=clb*0.52
    #vgm=clv*2.48
    #igm=cli*1.00

    #then convolve with moffat profile...

    bgm=scipy.signal.fftconvolve(bgm,psfb,mode='same')
    vgm=scipy.signal.fftconvolve(vgm,psfv,mode='same')
    igm=scipy.signal.fftconvolve(igm,psfi,mode='same')

    #save_to_jpeg(igm*5,vgm*5,bgm*5,2048,skyi,sigi,skyv,sigv,skyb,sigb,"slides_30.jpg")

    #arbitrary scaling? why? no idea...
    bgm=bgm*5.
    vgm=vgm*5.
    igm=igm*5.

    #rebin to desired pixel scale
    bgm=congrid.congrid(bgm,[1174,1174])
    vgm=congrid.congrid(vgm,[1174,1174])
    igm=congrid.congrid(igm,[1174,1174])

    save_to_jpeg(igm,vgm,bgm,1174,skyi,sigi,skyv,sigv,skyb,sigb,"slides_3.jpg")

    #NOW, add noise to each image to get desired noise scaling
    #scaling comes from measureing real SPT image
    nb=4.59
    nv=15.02
    ni=12.3

    bgm=bgm+np.random.normal(0,1,(1174L,1174L))*nb*0.95
    vgm=vgm+np.random.normal(0,1,(1174L,1174L))*nv*0.95
    igm=igm+np.random.normal(0,1,(1174L,1174L))*ni*0.95

    save_to_jpeg(igm,vgm,bgm,1174,skyi,sigi,skyv,sigv,skyb,sigb,"slides_4.jpg")

    #adds 'real' star images
    x1=200
    y1=1000

    bgm[x1-110:x1+110,y1-110:y1+110]=star1b*war1+nar1*bgm[x1-110:x1+110,y1-110:y1+110]
    vgm[x1-110:x1+110,y1-110:y1+110]=star1v*war1+nar1*vgm[x1-110:x1+110,y1-110:y1+110]
    igm[x1-110:x1+110,y1-110:y1+110]=star1i*war1+nar1*igm[x1-110:x1+110,y1-110:y1+110]

    x1=250
    y1=550

    bgm[x1-45:x1+45,y1-45:y1+45]=star2b*war2+nar2*bgm[x1-45:x1+45,y1-45:y1+45]
    vgm[x1-45:x1+45,y1-45:y1+45]=star2v*war2+nar2*vgm[x1-45:x1+45,y1-45:y1+45]
    igm[x1-45:x1+45,y1-45:y1+45]=star2i*war2+nar2*igm[x1-45:x1+45,y1-45:y1+45]

    x1=900
    y1=400

    igm[x1-43:x1+43,y1-43:y1+43]=star3i*war3+nar3*igm[x1-43:x1+43,y1-43:y1+43]
    vgm[x1-43:x1+43,y1-43:y1+43]=star3v*war3+nar3*vgm[x1-43:x1+43,y1-43:y1+43]
    bgm[x1-43:x1+43,y1-43:y1+43]=star3b*war3+nar3*bgm[x1-43:x1+43,y1-43:y1+43]

    ##add a little more noise to 'smooth' real stars into background
    #iim=iim+np.random.normal(0.0,1.0,(1174L,1174L))*ni/4.0
    #vim=vim+np.random.normal(0.0,1.0,(1174L,1174L))*nv/4.0
    #bim=bim+np.random.normal(0.0,1.0,(1174L,1174L))*nb/4.0

    save_to_jpeg(igm,vgm,bgm,1174,skyi,sigi,skyv,sigv,skyb,sigb,"slides_5.jpg")

    igm=igm+np.random.normal(0.0,1.0,(1174L,1174L))*ni/2.0
    vgm=vgm+np.random.normal(0.0,1.0,(1174L,1174L))*nv/2.0
    bgm=bgm+np.random.normal(0.0,1.0,(1174L,1174L))*nb/2.0

    save_to_jpeg(igm,vgm,bgm,1174,skyi,sigi,skyv,sigv,skyb,sigb,"slides_6.jpg")

    return igm,vgm,bgm

def cut_out_tails(im_in,dcut,ucut):

    im_in=im_in+2
    im_in[im_in<dcut] = dcut
    im_in[im_in>ucut] = ucut

    return im_in

def lensed_mapping_images_to_jpeg(iim,vim,bim,tag):

    skyi,sigi = plotim(iim)
    skyv,sigv = plotim(vim)
    skyb,sigb = plotim(bim)

    print skyi,sigi,skyv,sigv,skyb,sigb

    iim=(iim-skyi)/sigi
    vim=(vim-skyv)/sigv
    bim=(bim-skyb)/sigb

    acut=0.3
    icut=50.0*acut
    vcut=60.0*acut
    bcut=30.0*acut

    iim = cut_out_tails(iim,0,icut)
    vim = cut_out_tails(vim,0,vcut)
    bim = cut_out_tails(bim,0,bcut)

    iim=(iim*255/icut).astype(int)
    vim=(vim*255/vcut).astype(int)
    bim=(bim*255/bcut).astype(int)

    im=np.zeros((3,1174,1174),dtype=np.uint8)
    im[0,:,:]=iim
    im[1,:,:]=vim
    im[2,:,:]=bim

    scipy.misc.imsave('test_'+tag+'_'+'.jpg', im)

    return 0

def spt_mapping_images_to_jpeg():

    spti=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.i.fits')
    sptv=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.r.fits')
    sptb=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.g.fits')

    #sptb=sptb[5066-587:5066+587,5262-587:5262+587]
    ##sptv=sptv[5441-587:5441+587,5372-587:5372+587]
    #sptv=sptv[5066-587:5066+587,5262-587:5262+587]
    #spti=spti[5066-587:5066+587,5262-587:5262+587]

    spti=spti[5262-587:5262+587,5066-587:5066+587]
    sptv=sptv[5372-587:5372+587,5441-587:5441+587]
    sptb=sptb[5262-587:5262+587,5066-587:5066+587]

    skypi,sigpi=plotim(spti)
    skypv,sigpv=plotim(sptv)
    skypb,sigpb=plotim(sptb)

    print skypi,sigpi,skypv,sigpv,skypb,sigpb

    sptii=spti-skypi
    sptii=sptii/sigpi
    sptiv=sptv-skypv
    sptiv=sptiv/sigpv
    sptib=sptb-skypb
    sptib=sptib/sigpb

    acut=0.3
    icut=50.0*acut
    vcut=60.0*acut
    bcut=30.0*acut

    sptii = cut_out_tails(sptii,0,icut)
    sptiv = cut_out_tails(sptiv,0,vcut)
    sptib = cut_out_tails(sptib,0,bcut)

    sptii=(sptii*255/icut).astype(int)
    sptiv=(sptiv*255/vcut).astype(int)
    sptib=(sptib*255/bcut).astype(int)

    im=np.zeros((3,1174,1174),dtype=np.uint8)

    im[0,:,:]=sptii
    im[1,:,:]=sptiv
    im[2,:,:]=sptib

    #figure()
    #contour(sptii)
    #figure()
    #contour(sptiv)
    #figure()
    #contour(sptib)
    #show()

    im = np.uint8(im)

    scipy.misc.imsave('spt_comp.jpg', im)

    return 0

def save_to_jpeg(spti,sptv,sptb,nc,skypi,sigpi,skypv,sigpv,skypb,sigpb,filename):
    nx,ny = np.shape(spti)

    spti=spti[nx/2-nc/2:nx/2+nc/2,ny/2-nc/2:ny/2+nc/2]
    sptv=sptv[nx/2-nc/2:nx/2+nc/2,ny/2-nc/2:ny/2+nc/2]
    sptb=sptb[nx/2-nc/2:nx/2+nc/2,ny/2-nc/2:ny/2+nc/2]

    sptii=spti-skypi
    sptii=sptii/sigpi
    sptiv=sptv-skypv
    sptiv=sptiv/sigpv
    sptib=sptb-skypb
    sptib=sptib/sigpb

    acut=0.3
    icut=50.0*acut
    vcut=60.0*acut
    bcut=30.0*acut

    sptii = cut_out_tails(sptii,0,icut)
    sptiv = cut_out_tails(sptiv,0,vcut)
    sptib = cut_out_tails(sptib,0,bcut)

    sptii=(sptii*255/icut).astype(int)
    sptiv=(sptiv*255/vcut).astype(int)
    sptib=(sptib*255/bcut).astype(int)

    im=np.zeros((3,nc,nc),dtype=np.uint8)

    im[0,:,:]=sptii
    im[1,:,:]=sptiv
    im[2,:,:]=sptib

    im = np.uint8(im)

    scipy.misc.imsave(filename, im)

    return 0

def main():
    psf1,psf2,psfi,psfv,psfb = read_PSF_kernels()
    iim,vim,bim = read_lensed_images()

    iim = np.flipud(iim)*1.0
    vim = np.flipud(vim)*1.0
    bim = np.flipud(bim)*1.0

    #save_to_jpeg(iim,vim,bim,2048,"lensed_only.jpg")

    iom,vom,bom = read_los_images()

    #iim = (iim+iom)*1.0
    #vim = (vim+vom)*1.0
    #bim = (bim+bom)*1.0


    cli,clv,clb = read_mem_gals()

    star1i,star1v,star1b,war1,nar1,star2i,star2v,star2b,war2,nar2,star3i,star3v,star3b,war3,nar3 = grab_individual_bright_stars()
    #print star1i,war1,nar1
    igm,vgm,bgm = step_through_images(iim,vim,bim,iom*0.5,vom*0.5,bom*0.5,cli,clv,clb,star1i,star1v,star1b,star2i,star2v,star2b,star3i,star3v,star3b,war1,war2,war3,psfi,psfv,psfb)

    lensed_mapping_images_to_jpeg(igm,vgm,bgm,'test')

    #----------------------------------------------------------
    #spti=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.i.fits')
    #sptv=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.r.fits')
    #sptb=pyfits.getdata('swarp.SPT-CLJ0307-5042.2x2.g.fits')
    #
    #spt_mapping_images_to_jpeg()


    return 0

if __name__ == '__main__':
    main()
