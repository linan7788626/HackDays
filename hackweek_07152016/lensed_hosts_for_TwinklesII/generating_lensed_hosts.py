import numpy as np
import astropy.io.fits as pyfits
import om10_lensing_equations as ole
import subprocess as sp

apr = 206269.43 # arcseconds per rad

file_lagn_cat = "./data/sprinkled_agn_230.txt"
file_lens_cat = "./data/twinkles_lenses_v2.fits"
file_lgal_cat = "./data/sprinkled_lens_galaxies_230.txt"
file_host_cat = "./data/sprinkled_hosts_230.txt"

def create_inputs_for_ray_tracing_agnid(agnid):

    agn_id, agn_ra, agn_dec, agn_mag, agn_sed, agn_redshift, agn_twinkles_id, agn_img_num, agn_lens_galaxy_uID = np.loadtxt(
        file_lagn_cat, dtype="str", comments='#', delimiter=',',
        converters=None, skiprows=1, usecols=(1,2,3,4,5,6,17,18,19),
        unpack=True, ndmin=0)

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
    hdulist = pyfits.open(file_lens_cat)
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

#----------------------------------------------------------------------------
    lgal_id, lgal_ra, lgal_dec = np.loadtxt(file_lgal_cat, dtype="string",
                                            comments='#', delimiter=',',
                                            converters=None, skiprows=1,
                                            usecols=(1,2,3), unpack=True,
                                            ndmin=0)

    lgal_id = lgal_id.astype('int')
    idx_lgal = lgal_id == lens_galaxy_uID

    lens_ra = lgal_ra[idx_lgal][0]
    lens_dec = lgal_dec[idx_lgal][0]
#----------------------------------------------------------------------------
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
                'lens_ra'    : lens_ra,
                'lens_dec'   : lens_dec,
                'twinklesid' : twinklesid,
                'lensid'     : lid}

    # TTYPE8  = 'GAMMA'   TTYPE9  = 'PHIG    '
#----------------------------------------------------------------------------
# Parameters of the sources
#
    host_uid, host_ra, host_dec, host_mag, host_sed, host_redshift, \
    host_spatialmodel, host_major_axis, host_minor_axis, host_positionAngle, \
    host_sindex, host_twinkles_system, host_twinkles_img_num, \
    host_lens_galaxy_uid = np.loadtxt(file_host_cat, dtype="str", \
                                      comments='#', delimiter=',', \
                                      converters=None, skiprows=1, \
                                      usecols=(1,2,3,4,5,6,12,13,14,15,16,23,24,25), \
                                      unpack=True, ndmin=0)

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
# 2D Sersic Profile, no scale factor
#
def sersic_2d_hist(xi1,xi2,xc1,xc2,Ieff,Reff_arc,ql,pha,ndex=4.0):
    bn = 2.0*ndex-1/3.0+0.009876/ndex
    phirad = np.deg2rad(pha)
    xi1new = (xi1-xc1)*np.cos(phirad)+(xi2-xc2)*np.sin(phirad)
    xi2new = (xi2-xc2)*np.cos(phirad)-(xi1-xc1)*np.sin(phirad)
    R_scale = np.sqrt((xi1new**2)*ql+(xi2new**2)/ql)/Reff_arc
    res = np.exp(-bn*((R_scale)**(1.0/ndex)-1.0))
    return res
#--------------------------------------------------------------------
# Generate lensed Sersic Profile
#
def lensed_sersic(xi1,xi2,source_cat,lens_cat):

    xlc1 = lens_cat[0]          # x position of the lens, arcseconds
    xlc2 = lens_cat[1]          # y position of the lens, arcseconds
    rlc  = 0.0                  # core size of Non-singular Isothermal Ellipsoid
    # rle  = re_sv(lens_cat[2],lens_cat[5],source_cat[7])       #Einstein radius of lens, arcseconds.
    rle = ole.re_sv(lens_cat[2],lens_cat[5],source_cat[7])

    ql   = lens_cat[3]          # axis ratio b/a
    phl = lens_cat[4]           # orientation, degree
    ext_shears = lens_cat[6]
    ext_angle = lens_cat[7]
    #----------------------------------------------------------------------
    le = ole.lambda_e_tot(1.0-ql)

    ai1, ai2 = ole.alphas_sie(xlc1, xlc2, phl, ql, rle, le, ext_shears, ext_angle, rlc, xi1, xi2)

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

    g_source = sersic_2d_hist(xi1,xi2,ysc1,ysc2,1.0,Reff_arc,qs,phs,ndex)
    g_limage = sersic_2d_hist(yi1,yi2,ysc1,ysc2,1.0,Reff_arc,qs,phs,ndex)

    mag_lensed = mag_tot - 2.5*np.log(np.sum(g_limage)/np.sum(g_source))

    return mag_lensed, g_limage


def generate_lensed_host(agnID):

    lensP, srcsP = create_inputs_for_ray_tracing_agnid(agnID)

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

        lens_ra = lensP['lens_ra']
        lens_dec = lensP['lens_dec']

        file_limg = "./output_fits/"+str(agnID)+"_"+str(i)+"_"+str(lens_ra)+"_"+str(lens_dec)+"_"+str(lensed_mag)+"_"+str(srcsP[i]['sed_src'].split("/")[0])+"_"+str(srcsP[i]['sed_src'].split("/")[1])+"_"+str(srcsP[i]['zs'])+"_"+str(dsx)+".fits"

        pyfits.writeto(file_limg, lensed_image.astype("float32"), overwrite=True)

        cmd = "bzip2 -f " + file_limg
        sp.call(cmd, shell=True)

    return 0


if __name__ == '__main__':
    agn_id, agn_img_num = np.loadtxt(file_lagn_cat, dtype="str", comments='#', delimiter=',', converters=None, skiprows=1, usecols=(1,18), unpack=True, ndmin=0)

    agn_id = agn_id.astype('int')
    agn_img_num = agn_img_num.astype('int')

    img_num = 0
    idx = agn_img_num==img_num
    AGN_ID = agn_id[idx]

    for i in AGN_ID:
        print i
        generate_lensed_host(i)
