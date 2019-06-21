import numpy as np
from astropy.cosmology import Planck13 as p13
from astropy import constants as const
import astropy.io.fits as pyfits
import scipy.special as spf
# import subprocess as sp

apr = 206269.43 # arcseconds per rad

def create_inputs_for_ray_tracing_agnid(agnid):
    agn_id, agn_ra, agn_dec, agn_mag, agn_sed, agn_redshift, agn_catsim_id, agn_twinkles_id, agn_img_num = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230_catids.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)

    agn_id = agn_id.astype('int')
    idx_agn = agn_id == agnid

    agn_catsim_id = agn_catsim_id.astype('int')
    agn_twinkles_id = agn_twinkles_id.astype('int')
    agn_img_num = agn_img_num.astype('int')
#----------------------------------------------------------------------------
# prefix,uniqueId,raPhoSim,decPhoSim,phoSimMagNorm,sedFilepath,redshift,shear1,shear2,kappa,raOffset,decOffset,spatialmodel,majorAxis,minorAxis,positionAngle,sindex,internalExtinctionModel,internalAv,internalRv,galacticExtinctionModel,galacticAv,galacticRv,twinkles_system,twinkles_img_num,lens_galaxy_uID
    host_catsim_id, host_ra, host_dec, host_mag, host_sed, host_redshift, host_spatialmodel, host_major_axis, host_minor_axis, host_positionAngle, host_sindex, host_id = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_hosts_230.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(1,2,3,4,5,6,12,13,14,15,16,26), unpack=True, ndmin=0)

    host_id = host_id.astype('int')
    host_catsim_id = host_catsim_id.astype('int')
    host_mag = host_mag.astype('double')
    host_redshift = host_redshift.astype('double')
    host_major_axis = host_major_axis.astype('double')
    host_minor_axis = host_minor_axis.astype('double')
    host_positionAngle = host_positionAngle.astype('double')
    host_sindex = host_sindex.astype('int')
    host_catsim_id = host_catsim_id.astype('int')
#----------------------------------------------------------------------------
# Parameters of the lens, by matching twinklesid
#
    twinklesid = agn_twinkles_id[idx_agn][0]
    hdulist = pyfits.open('./twinkles_DESC_SLAC/twinkles_lenses_v2.fits')
    # lid = hdulist[1].data['LENSID'][twinklesid]
    xl1 = 0.0
    xl2 = 0.0
    vd = hdulist[1].data['VELDISP'][twinklesid]   # needed from OM10
    zd = hdulist[1].data['ZLENS'][twinklesid]
    ql  = 1.0 - hdulist[1].data['ELLIP'][twinklesid]
    phi= hdulist[1].data['PHIE'][twinklesid]

    ys1 = hdulist[1].data['XSRC'][twinklesid]
    ys2 = hdulist[1].data['YSRC'][twinklesid]

    ext_shr = hdulist[1].data['GAMMA'][twinklesid]
    ext_phi = hdulist[1].data['PHIG'][twinklesid]

    lens_cat = {'xl1' : xl1, 'xl2' : xl2, 'ql' : ql, 'vd' : vd,
                'phl' : phi, 'gamma' : ext_shr, 'phg':ext_phi,
                'zl' : zd, 'twinklesid' : twinklesid}

    # TTYPE8  = 'GAMMA'   TTYPE9  = 'PHIG    '
#----------------------------------------------------------------------------
# Parameters of the sources
#
    idx2 = agn_img_num==0
    idx = idx_agn&idx2

    catsim_id = agn_catsim_id[idx]
    idx_host = host_catsim_id == catsim_id
    hostid = host_id[idx_host]
    srcs_cats = []
    for i in xrange(len(hostid)):
        srcs_cat = {}

        idx_host_components =  host_id == hostid[i]

        mag_src = host_mag[idx_host_components][0]
        Reff_src = np.sqrt(host_major_axis[idx_host_components][0]
                           *host_minor_axis[idx_host_components][0])
        qs = host_minor_axis[idx_host_components][0]/host_major_axis[idx_host_components][0]
        phs = host_positionAngle[idx_host_components][0]
        ns = host_sindex[idx_host_components][0]
        zs = host_redshift[idx_host_components][0]
        sed_src = host_sed[idx_host_components][0]

        srcs_cat = {'ys1' : ys1, 'ys2' : ys2, 'mag_src' : mag_src,
                    'Reff_src' : Reff_src, 'qs' : qs, 'phs' : phs,
                    'ns' : ns, 'zs' : zs, 'sed_src' : sed_src,
                    'agnid' : agnid, 'componentsid' : i}

        srcs_cats.append(srcs_cat)

        # lensed_images, lensed_mag, lensed_sed = ray_tracing(lens_cat, srcs_cat)
        # the output should be an image with TwinkleID + SED + Mag + zl +

    return lens_cat, srcs_cats
#--------------------------------------------------------------------
# Construct regular grids
#
def make_r_coor(nc,dsx):

    bsz = nc*dsx
    x1 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0
    x2 = np.linspace(0,bsz-dsx,nc)-bsz/2.0+dsx/2.0

    x1,x2 = np.meshgrid(x1,x2)
    return x1,x2
#--------------------------------------------------------------------
# Calculate deflection angles
#
def lens_equation_sie(x1,x2,xc1,xc2,q,rc,re,pha):

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
# Rotate regular grids
#
def xy_rotate(x, y, xcen, ycen, phi):
    phirad = np.deg2rad(phi)
    xnew = (x-xcen)*np.cos(phirad)+(y-ycen)*np.sin(phirad)
    ynew = (y-ycen)*np.cos(phirad)-(x-xcen)*np.sin(phirad)
    return (xnew,ynew)
#--------------------------------------------------------------------
# Calculate Einstein Radius according to Velocity Dispersion
#
def re_sv(sv,z1,z2):
    Da_s = p13.angular_diameter_distance(z2).value
    Da_ls = p13.angular_diameter_distance_z1z2(z1,z2).value
    res = 4.0*np.pi*(sv**2.0/(const.c.value/1e3)**2.0)*Da_ls/Da_s*apr
    return res
#--------------------------------------------------------------------
# 2D Sersic Profile
#
def sersic_2d(xi1,xi2,xc1,xc2,Ieff,Reff_arc,ql,pha,ndex=4.0):
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    (xi1new,xi2new) = xy_rotate(xi1, xi2, xc1, xc2, pha)
    R_scale = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff_arc
    res = Ieff*np.exp(-bn*((R_scale)**(1.0/ndex)-1.0))
    return res
#--------------------------------------------------------------------
# Convert total magnitude to effective flux
#
def sersic_mag_tot_to_Ieff(mag_tot,Reff,ndex,z,Rcut=100): # g band
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    xtmp = bn*(Rcut)**(1.0/ndex)
    ktmp = 2.0*np.pi*ndex*(np.e**bn/(bn**(2.0*ndex)))*spf.gammainc(2.0*ndex,xtmp)*spf.gamma(2.0*ndex)
    Dl_s = p13.luminosity_distance(z).value*1e6

    Ieff = 10.0**((5.12-2.5*np.log10(ktmp)-5.0*np.log10(Reff)+5.0*np.log10(Dl_s/10)-mag_tot)*0.4) # L_sun/pc^2
    return Ieff
#--------------------------------------------------------------------
# Generate lensed Sersic Profile
#
def lensed_sersic(xi1,xi2,source_cat,lens_cat):

    xlc1 = lens_cat[0]          # x position of the lens, arcseconds
    xlc2 = lens_cat[1]          # y position of the lens, arcseconds
    rlc  = 0.0                  # core size of Non-singular Isothermal Ellipsoid
    rle  = re_sv(lens_cat[2],lens_cat[5],source_cat[7])       #Einstein radius of lens, arcseconds.
    ql   = lens_cat[3]          # axis ratio b/a
    phl = lens_cat[4]           # orientation, degree
    #----------------------------------------------------------------------
    ai1,ai2,mua = lens_equation_sie(xi1,xi2,xlc1,xlc2,ql,rlc,rle,phl)

    yi1 = xi1-ai1
    yi2 = xi2-ai2
    #----------------------------------------------------------------------
    ysc1 = source_cat[0]        # x position of the source, arcseconds
    ysc2 = source_cat[1]        # y position of the source, arcseconds
    #----------------------------------------------------------------------
    mag_tot = source_cat[2]     # total magnitude of the source
    Reff_arc = source_cat[3]    # Effective Radius of the source, arcseconds
    qs   = source_cat[4]        # axis ratio of the source, b/a
    phs  = source_cat[5]        # orientation of the source, degree
    ndex = source_cat[6]        # index of the source
    zs = source_cat[7]          # redshift of the source

    Da_s = p13.angular_diameter_distance(zs).value
    Reff = Reff_arc*Da_s/apr*1e6 # pc
    Ieff = sersic_mag_tot_to_Ieff(mag_tot,Reff,ndex,zs) #L_sun/pc^2
    g_limage = sersic_2d(yi1,yi2,ysc1,ysc2,Ieff,Reff_arc,qs,phs,ndex)
    g_source = sersic_2d(xi1,xi2,ysc1,ysc2,Ieff,Reff_arc,qs,phs,ndex)

    mag_lensed = mag_tot - 2.5*np.log(np.sum(g_limage)/np.sum(g_source))

    return mag_lensed, g_limage


def generate_lensed_host(agnID):
    lensP, srcsP = create_inputs_for_ray_tracing_agnid(agnID)

    for i in xrange(len(srcsP)):
        #-----------------------------------------------------------------------
        # parameters of the mass model of the lens
        #
        lens_cats = [lensP['xl1'], lensP['xl2'], lensP['vd'], lensP['ql'], lensP['phl'], lensP['zl']]
        #-----------------------------------------------------------------------
        # parameters of the host galaxy of lensed point source
        #
        srcs_cats = [srcsP[i]['ys1'], srcsP[i]['ys2'], srcsP[i]['mag_src'], srcsP[i]['Reff_src'], srcsP[i]['qs'], srcsP[i]['phs'], srcsP[i]['ns'], srcsP[i]['zs']]
        #-----------------------------------------------------------------------
        # Calculate the lensed images
        #
        dsx = 0.01  # pixel size per side, arcseconds
        nnn = 1000  # number of pixels per side
        xi1, xi2 = make_r_coor(nnn,dsx)
        #-----------------------------------------------------------------------
        # Calculate the lensed image and its total magnitude
        #
        lensed_mag, lensed_image = lensed_sersic(xi1,xi2,srcs_cats,lens_cats)
        #-----------------------------------------------------------------------
        # Visualize the lensed imagesr. If we normalize the image, the 2D matrix
        # becomes the normalized histogram of the light distribution.
        #
        agn_id_tmp, agn_lensgal_id_tmp = np.loadtxt("./twinkles_DESC_SLAC/filter_0_sprinkled_AGNs.txt", dtype="str", comments='#', delimiter=' ', converters=None, skiprows=1, usecols=(0,13), unpack=True, ndmin=0)

        agn_id_tmp = agn_id_tmp.astype('int')
        idx_agn_tmp = agn_id_tmp == agnID

        agn_lensgal_id_tmp = agn_lensgal_id_tmp.astype('int')
#----------------------------------------------------------------------------
        lgal_id, lgal_ra, lgal_dec, lgal_mag_norm, lgal_redshift, lgal_major_axis, lgal_minor_axis = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_lens_galaxies_230.txt", dtype="string", comments='#', delimiter=',', converters=None, skiprows=1, usecols=None, unpack=True, ndmin=0)
        lgal_id = lgal_id.astype('int')

        idx_lgal = lgal_id == agn_lensgal_id_tmp[idx_agn_tmp]

        lens_ra = lgal_ra[idx_lgal][0]
        lens_dec = lgal_dec[idx_lgal][0]

        file_limg = "./outputs_fits/"+str(agnID)+"_"+str(i)+"_"+str(lens_ra)+"_"+str(lens_dec)+"_"+str(lensed_mag)+"_"+str(srcsP[i]['sed_src'].split("/")[0])+"_"+str(srcsP[i]['sed_src'].split("/")[1])+"_"+str(srcsP[i]['zs'])+"_"+str(dsx)+".fits"
        pyfits.writeto(file_limg, lensed_image.astype("float32"), overwrite=True)

        # cmd = "bzip2 -f " + file_limg
        # sp.call(cmd, shell=True)

    return 0

if __name__ == '__main__':
    AGN_ID = np.loadtxt("./twinkles_DESC_SLAC/filter_0_sprinkled_AGNs.txt", dtype="int", comments='#', delimiter=' ', converters=None, skiprows=1, usecols=(0,), unpack=True, ndmin=0)

    for i in AGN_ID:
        print i
        generate_lensed_host(i)

