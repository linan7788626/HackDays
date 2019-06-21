import numpy as np
from astropy.cosmology import Planck13 as p13
from astropy import constants as const
import astropy.io.fits as pyfits
import scipy.special as spf
import om10_lensing_equations as ole
import lambda_e as lef
import pylab as pl

import om10
db = om10.DB(catalog="/Users/uranus/GitHub/OM10/data/qso_mock.fits")

# import subprocess as sp

apr = 206269.43 # arcseconds per rad

def create_inputs_for_ray_tracing_agnid(agnid):
    # prefix,uniqueId,raPhoSim,decPhoSim,phoSimMagNorm,sedFilepath,redshift,shear1,shear2,kappa,raOffset,decOffset,spatialmodel,internalExtinctionModel,galacticExtinctionModel,galacticAv,galacticRv,twinkles_system,twinkles_img_num,lens_galaxy_uID
    agn_id, agn_ra, agn_dec, agn_mag, agn_sed, agn_redshift, agn_twinkles_id, agn_img_num, agn_lens_galaxy_uID = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_agn_230.txt",
               dtype="str", comments='#', delimiter=',', converters=None,
               skiprows=1, usecols=(1,2,3,4,5,6,17,18,19), unpack=True,
               ndmin=0)

    agn_id = agn_id.astype('int')
    idx_agn = agn_id == agnid

    agn_twinkles_id = agn_twinkles_id.astype('int')
    agn_img_num = agn_img_num.astype('int')
    agn_lens_galaxy_uID = agn_lens_galaxy_uID.astype('int')
    agn_ra = agn_ra.astype('double')
    agn_dec = agn_dec.astype('double')

    lens_galaxy_uID = agn_lens_galaxy_uID[idx_agn]
#----------------------------------------------------------------------------
# Parameters of the lens, by matching twinklesid
#
    twinklesid = agn_twinkles_id[idx_agn][0]
    hdulist = pyfits.open('./twinkles_DESC_SLAC/twinkles_lenses_v2.fits')
    # print hdulist[1].header   # needed from OM10

    idx = hdulist[1].data['twinklesId'] == twinklesid
    lid = hdulist[1].data['LENSID'][idx][0]
    xl1 = 0.0
    xl2 = 0.0
    vd = hdulist[1].data['VELDISP'][idx][0]   # needed from OM10
    zd = hdulist[1].data['ZLENS'][idx][0]
    ql  = 1.0 - hdulist[1].data['ELLIP'][idx][0]
    phi= hdulist[1].data['PHIE'][idx][0]

    ys1 = hdulist[1].data['XSRC'][idx][0]
    ys2 = hdulist[1].data['YSRC'][idx][0]

    ext_shr = hdulist[1].data['GAMMA'][idx][0]
    ext_phi = hdulist[1].data['PHIG'][idx][0]

    ximg = hdulist[1].data['XIMG'][idx][0]
    yimg = hdulist[1].data['YIMG'][idx][0]

    lens_cat = {'xl1'        : xl1,
                'xl2'        : xl2,
                'ql'         : ql,
                'vd'         : vd,
                'phl'        : phi,
                'gamma'      : ext_shr,
                'phg'        : ext_phi,
                'zl'         : zd,
                'ximg'       : ximg,
                'yimg'       : yimg,
                'twinklesid' : twinklesid,
                'lensid'     : lid}

    # TTYPE8  = 'GAMMA'   TTYPE9  = 'PHIG    '
#----------------------------------------------------------------------------
# Parameters of the sources
#

# prefix,uniqueId,raPhoSim,decPhoSim,phoSimMagNorm,sedFilepath,redshift,shear1,shear2,kappa,raOffset,decOffset,spatialmodel,majorAxis,minorAxis,positionAngle,sindex,internalExtinctionModel,internalAv,internalRv,galacticExtinctionModel,galacticAv,galacticRv,twinkles_system,twinkles_img_num,lens_galaxy_uID
    host_uid, host_ra, host_dec, host_mag, host_sed, host_redshift, host_spatialmodel, host_major_axis, host_minor_axis, host_positionAngle, host_sindex, host_twinkles_system, host_twinkles_img_num, host_lens_galaxy_uid = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_hosts_230.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(1,2,3,4,5,6,12,13,14,15,16,23,24,25), unpack=True, ndmin=0)

    host_uid = host_uid.astype('int')
    host_mag = host_mag.astype('double')
    host_ra = host_ra.astype('double')
    host_dec = host_dec.astype('double')
    host_redshift = host_redshift.astype('double')
    host_major_axis = host_major_axis.astype('double')*apr
    host_minor_axis = host_minor_axis.astype('double')*apr
    host_positionAngle = host_positionAngle.astype('double')
    host_sindex = host_sindex.astype('int')
    host_twinkles_system = host_twinkles_system.astype('int')
    host_twinkles_img_num = host_twinkles_img_num.astype('int')
    host_lens_galaxy_uid = host_lens_galaxy_uid.astype('int')

    idx_host1 = host_lens_galaxy_uid == lens_galaxy_uID
    idx_host2 = host_twinkles_img_num == 0
    idx_host = idx_host1&idx_host2

    srcs_cats = []
    for i in xrange(len(host_mag[idx_host])):
        srcs_cat = {}

        mag_src = host_mag[idx_host][i]
        Reff_src = np.sqrt(host_major_axis[idx_host][i]
                           *host_minor_axis[idx_host][i])
        qs = host_minor_axis[idx_host][i]/host_major_axis[idx_host][i]
        phs = host_positionAngle[idx_host][i]
        ns = host_sindex[idx_host][i]
        zs = host_redshift[idx_host][i]
        # zs = hdulist[1].data['ZSRC'][idx][0]
        sed_src = host_sed[idx_host][i]

        srcs_cat = {'ys1'          : ys1,
                    'ys2'          : ys2,
                    'mag_src'      : mag_src,
                    'Reff_src'     : Reff_src,
                    'qs'           : qs,
                    'phs'          : phs,
                    'ns'           : ns,
                    'zs'           : zs,
                    'sed_src'      : sed_src,
                    'agnid'        : agnid,
                    'ra'           : host_ra[idx_host][i],
                    'dec'          : host_dec[idx_host][i],
                    'lensed_agn_ra': agn_ra[idx_agn][0],
                    'lensed_agn_dec': agn_dec[idx_agn][0],
                    'lensid'       : host_lens_galaxy_uid[idx_host][i],
                    'componentsid' : i}

        srcs_cats.append(srcs_cat)

    # lensed_images, lensed_mag, lensed_sed = ray_tracing(lens_cat, srcs_cat)
    # the output should be an image with TwinkleID + SED + Mag + zl +
#----------------------------------------------------------------------------
    #prefix,uniqueId,raPhoSim,decPhoSim,phoSimMagNorm,sedFilepath,redshift,shear1,shear2,kappa,raOffset,decOffset,spatialmodel,majorAxis,minorAxis,positionAngle,sindex,internalExtinctionModel,internalAv,internalRv,galacticExtinctionModel,galacticAv,galacticRv
    lens_gals_uid, lens_gals_ra, lens_gals_dec, lens_gals_mag, lens_gals_sed, lens_gals_redshift, lens_gals_spatialmodel, lens_gals_major_axis, lens_gals_minor_axis, lens_gals_positionAngle, lens_gals_sindex = np.loadtxt("./twinkles_DESC_SLAC/sprinkled_lens_galaxies_230.txt", dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(1,2,3,4,5,6,12,13,14,15,16), unpack=True, ndmin=0)

    lens_gals_uid = lens_gals_uid.astype('int')
    lens_gals_mag = lens_gals_mag.astype('double')
    lens_gals_ra = lens_gals_ra.astype('double')
    lens_gals_dec = lens_gals_dec.astype('double')
    lens_gals_redshift = lens_gals_redshift.astype('double')
    lens_gals_major_axis = lens_gals_major_axis.astype('double')
    lens_gals_minor_axis = lens_gals_minor_axis.astype('double')
    lens_gals_positionAngle = lens_gals_positionAngle.astype('double')
    lens_gals_sindex = lens_gals_sindex.astype('int')

    idx_lens_gals = lens_gals_uid == lens_galaxy_uID

    lens_gals_cat = {'ra'         : lens_gals_ra[idx_lens_gals][0],
                     'dec'        : lens_gals_dec[idx_lens_gals][0],
                     'mag'        : lens_gals_mag[idx_lens_gals][0],
                     'a'          : lens_gals_major_axis[idx_lens_gals][0],
                     'b'          : lens_gals_minor_axis[idx_lens_gals][0],
                     'phi'        : lens_gals_positionAngle[idx_lens_gals][0],
                     'smodel'     : lens_gals_spatialmodel[idx_lens_gals][0],
                     'sindex'     : lens_gals_sindex[idx_lens_gals][0],
                     'zl'         : lens_gals_redshift[idx_lens_gals][0],
                     'twinklesid' : twinklesid,
                     'lensid'     : lens_galaxy_uID[0]}
#--------------------------------------------------------------------
    return lens_gals_cat, lens_cat, srcs_cats
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
    # rle  = re_sv(lens_cat[2],lens_cat[5],source_cat[7])       #Einstein radius of lens, arcseconds.
    rle = ole.re_sv(lens_cat[2],lens_cat[5],source_cat[7])

    ql   = lens_cat[3]          # axis ratio b/a
    phl = lens_cat[4]           # orientation, degree
    ext_shears = lens_cat[6]
    ext_angle = lens_cat[7]
    #----------------------------------------------------------------------
    le = lef.lambda_e_tot(1.0-ql)

    ai1, ai2 = ole.alphas_sie(xlc1, xlc2, phl, ql, rle, le, ext_shears, ext_angle, 0.0, xi1, xi2)

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
    lgalP, lensP, srcsP = create_inputs_for_ray_tracing_agnid(agnID)

    for i in xrange(len(srcsP)):
        #-----------------------------------------------------------------------
        # parameters of the mass model of the lens
        #
        lens_cats = [lensP['xl1'], lensP['xl2'], lensP['vd'], lensP['ql'], lensP['phl'], lensP['zl'], lensP['gamma'], lensP['phg']]
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
#----------------------------------------------------------------------------
        print lensP['ximg']
        print lensP['yimg']
        print (srcsP[i]['lensed_agn_ra']-lgalP['ra'])/np.pi*180/3600.0*apr, (srcsP[i]['lensed_agn_dec']-lgalP['dec'])/np.pi*180/3600.0*apr

        # l_id = lensP['lensid']
        # print l_id
        # lens_tmp = db.get_lens(l_id)
        # om10.plot_lens(lens_tmp)

        pl.figure(figsize=(8,8))
        pl.contourf(xi1,xi2,lensed_image)
        pl.plot(lensP['ximg'][np.nonzero(lensP['ximg'])], \
                lensP['yimg'][np.nonzero(lensP['yimg'])], 'bx')

        # pl.savefig("./outputs_pngs/"+str(agnID)+"_"+str(i)+".png")

        # lens_ra = lgalP['ra']
        # lens_dec = lgalP['dec']

        # file_limg = "./outputs_fits/"+str(srcsP[i]['lensid'])+"_"+str(i)+"_"+str(lens_ra)+"_"+str(lens_dec)+"_"+str(lensed_mag)+"_"+str(srcsP[i]['sed_src'].split("/")[0])+"_"+str(srcsP[i]['sed_src'].split("/")[1])+"_"+str(srcsP[i]['zs'])+"_"+str(dsx)+".fits"
        # pyfits.writeto(file_limg, lensed_image.astype("float32"), overwrite=True)

        # cmd = "bzip2 -f " + file_limg
        # sp.call(cmd, shell=True)

    return 0

if __name__ == '__main__':
    AGN_ID = np.loadtxt("./twinkles_DESC_SLAC/filter_0_sprinkled_AGNs.txt", dtype="int", comments='#', delimiter=' ', converters=None, skiprows=1, usecols=(0,), unpack=True, ndmin=0)

    for i in AGN_ID[:10]:
        print i
        generate_lensed_host(i)

    pl.show()
