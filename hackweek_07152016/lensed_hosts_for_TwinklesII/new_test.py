import numpy as np
import astropy.io.fits as pyfits
import om10_lensing_equations as ole
import subprocess as sp
import pickle

apr = 206269.43 # arcseconds per rad

file_cat = "./data/glsn_protodc2.pkl"
# file_lagn_cat = "./data/sprinkled_agn_230.txt"
# file_lens_cat = "./data/twinkles_lenses_v2.fits"
# file_lgal_cat = "./data/sprinkled_lens_galaxies_230.txt"
# file_host_cat = "./data/sprinkled_hosts_230.txt"

def create_inputs_for_ray_tracing_agnid(index, input_cats):

#----------------------------------------------------------------------------
# Parameters of the lens, by matching twinklesid
#

    sysid = index
    xl1 = 0.0
    xl2 = 0.0
    vd  = input_cats['sigma'][index]   # needed from OM10
    zd  = input_cats['zl'][index]
    ql  = 1.0 - input_cats['e'][index]
    phi = input_cats['theta_e'][index]

    yhs1 = input_cats['xhost'][index]
    yhs2 = input_cats['yhost'][index]

    yns1 = input_cats['snx'][index]
    yns2 = input_cats['sny'][index]

    ext_shr = input_cats['gamma'][index]
    ext_phi = input_cats['theta_gamma'][index]

#----------------------------------------------------------------------------
    lens_cat = {'xl1'        : xl1,
                'xl2'        : xl2,
                'ql'         : ql,
                'vd'         : vd,
                'phl'        : phi,
                'gamma'      : ext_shr,
                'phg'        : ext_phi,
                'zl'         : zd,
                'sysid'      : sysid}

    # TTYPE8  = 'GAMMA'   TTYPE9  = 'PHIG    '
#----------------------------------------------------------------------------
# Parameters of the sources
#

    srcs_cats = []
    num_components = 2
    for i in xrange(num_components):
        srcs_cat = {}

        mag_src = host_mag[idx_host][i]
        Reff_src = np.sqrt(host_major_axis[idx_host][i]
                           *host_minor_axis[idx_host][i])
        qs = host_minor_axis[idx_host][i]/host_major_axis[idx_host][i]
        phs = host_positionAngle[idx_host][i]
        ns = host_sindex[idx_host][i]
        zs = host_redshift[idx_host][i]
        sed_src = host_sed[idx_host][i]

        srcs_cat = {'yhs1'         : yhs1,
                    'yhs2'         : yhs2,
                    'yns1'         : yns1,
                    'yns2'         : yhs2,
                    'mag_src'      : mag_src,
                    'Reff_src'     : Reff_src,
                    'qs'           : qs,
                    'phs'          : phs,
                    'ns'           : ns,
                    'zs'           : zs,
                    'sed_src'      : sed_src,
                    'sysid'        : sysid,
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

    input_cat = pickle.load(open('./data/glsn_protodc2.pkl','rb'))
    # num_sims = len(input_cat['galaxy_id'])
    num_sims = 1

    for i in xrange(num_sims):
        print i
        generate_lensed_host(i, input_cat)
