import numpy as np
import ctypes as ct
#--------------------------------------------------------------------
icic = ct.CDLL("../improving_OM10_using_Illustris/lib_so_icic/libicic.so")
icic.inverse_cic.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_double),\
                            np.ctypeslib.ndpointer(dtype =  ct.c_double), \
                            np.ctypeslib.ndpointer(dtype =  ct.c_double), \
                            ct.c_double,ct.c_double,ct.c_double,ct.c_int,ct.c_int,ct.c_int,\
                            np.ctypeslib.ndpointer(dtype = ct.c_double)]
icic.inverse_cic.restype  = ct.c_void_p

def call_inverse_cic_single(img_in,yi1,yi2,dsi):
    yc1 = 0.0
    yc2 = 0.0
    ny1,ny2 = np.shape(img_in)
    nlimgs = len(yi1)

    img_in = np.array(img_in,dtype=ct.c_double)

    yi1 = np.array(yi1,dtype=ct.c_double)
    yi2 = np.array(yi2,dtype=ct.c_double)

    img_out = np.zeros((nlimgs))

    icic.inverse_cic(img_in,yi1,yi2,ct.c_double(yc1),ct.c_double(yc2),ct.c_double(dsi),ct.c_int(ny1),ct.c_int(ny2),ct.c_int(nlimgs),img_out)
    return img_out

rtf = ct.CDLL("../improving_OM10_using_Illustris/lib_so_omp_icic/libicic.so")
rtf.inverse_cic.argtypes = [np.ctypeslib.ndpointer(dtype =  ct.c_double),\
                            np.ctypeslib.ndpointer(dtype =  ct.c_double), \
                            np.ctypeslib.ndpointer(dtype =  ct.c_double), \
                            ct.c_double,ct.c_double,ct.c_double,ct.c_int,ct.c_int,ct.c_int,ct.c_int,\
                            np.ctypeslib.ndpointer(dtype = ct.c_double)]
rtf.inverse_cic.restype  = ct.c_void_p

def call_inverse_cic(img_in,yi1,yi2,dsi):
    yc1 = 0.0
    yc2 = 0.0
    ny1,ny2 = np.shape(img_in)
    nx1,nx2 = np.shape(yi1)

    img_in = np.array(img_in,dtype=ct.c_double)

    yi1 = np.array(yi1,dtype=ct.c_double)
    yi2 = np.array(yi2,dtype=ct.c_double)

    img_out = np.zeros((nx1,nx2))

    rtf.inverse_cic(img_in,yi1,yi2,ct.c_double(yc1),ct.c_double(yc2),ct.c_double(dsi),ct.c_int(ny1),ct.c_int(ny2),ct.c_int(nx1),ct.c_int(nx2),img_out)
    return img_out.reshape((nx1,nx2))

#--------------------------------------------------------------------
gls = ct.CDLL("../improving_OM10_using_Illustris/lib_so_cal_alphas/lib/libglsg.so")

gls.sdens_to_shears.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
								ct.c_int,ct.c_double,ct.c_double,ct.c_double,\
								np.ctypeslib.ndpointer(dtype = ct.c_double), \
								np.ctypeslib.ndpointer(dtype = ct.c_double)]
gls.sdens_to_alphas.restype  = ct.c_void_p

def call_cal_shears(sdens,Ncc,boxsize,zl,zs):

	sdens = np.array(sdens,dtype=ct.c_double)
	shear1 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)
	shear2 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)

	gls.sdens_to_shears(sdens,Ncc,boxsize,zl,zs,shear1,shear2)

	return shear1,shear2

gls.sdens_to_alphas.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
								ct.c_int,ct.c_double,ct.c_double,ct.c_double,\
								np.ctypeslib.ndpointer(dtype = ct.c_double), \
								np.ctypeslib.ndpointer(dtype = ct.c_double)]
gls.sdens_to_alphas.restype  = ct.c_void_p

def call_cal_alphas0(sdens,Ncc,boxsize,zl,zs):

	sdens = np.array(sdens,dtype=ct.c_double)
	alpha1 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)
	alpha2 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)

	gls.sdens_to_alphas(sdens,Ncc,boxsize,zl,zs,alpha1,alpha2)

	return alpha1,alpha2

gls.nkappa_to_alphas.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                 np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                 np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                 ct.c_int,ct.c_double]
gls.nkappa_to_alphas.restype  = ct.c_void_p

def call_cal_alphas1(kappa,Ncc,Dcell):

    kappa = kappa.astype(ct.c_double)
    alpha1 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)
    alpha2 = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)

    gls.nkappa_to_alphas(kappa,alpha1,alpha2,ct.c_int(Ncc),ct.c_double(Dcell))
    return alpha1,alpha2

gls.alphas_to_mu.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                 np.ctypeslib.ndpointer(dtype = ct.c_double), \
                                 ct.c_int,ct.c_double, \
                             np.ctypeslib.ndpointer(dtype = ct.c_double)]
gls.alphas_to_mu.restype  = ct.c_void_p

def call_alphas_to_mu(alpha1,alpha2,Ncc,Dcell):

    alpha1 = np.array(alpha1,dtype=ct.c_double)
    alpha2 = np.array(alpha2,dtype=ct.c_double)
    mu = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)

    gls.alphas_to_mu(alpha1,alpha2,ct.c_int(Ncc),ct.c_double(Dcell),mu)
    return mu

gls.sdens_to_phi.argtypes = [np.ctypeslib.ndpointer(dtype = ct.c_double), \
								ct.c_int,ct.c_double,ct.c_double,ct.c_double,\
								np.ctypeslib.ndpointer(dtype = ct.c_double)]
gls.sdens_to_phi.restype  = ct.c_void_p

def call_cal_phi(sdens,Ncc,boxsize,zl,zs):

	sdens = np.array(sdens,dtype=ct.c_double)
	phi = np.array(np.zeros((Ncc,Ncc)),dtype=ct.c_double)
	gls.sdens_to_phi(sdens,Ncc,boxsize,zl,zs,phi)

	return phi
#--------------------------------------------------------------------
sps = ct.CDLL("../improving_OM10_using_Illustris/lib_so_cal_sdens_sph/lib/libsphsdens.so")
sps.cal_sph_sdens_arrays.argtypes =[np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     ct.c_float,ct.c_long, \
                                     ct.c_float,ct.c_long,ct.c_long, \
                                     ct.c_float,ct.c_float,ct.c_float,ct.c_float, \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float)]
sps.cal_sph_sdens_arrays.restype  = ct.c_int

def call_sph_sdens_arrays(xi1,xi2,xi3,xic1,xic2,xic3,Bsz,Nc,Np,pmass):

    x1 = xi1.astype(ct.c_float)
    x2 = xi2.astype(ct.c_float)
    x3 = xi3.astype(ct.c_float)
    #x1 = np.array(xi1,dtype=ct.c_float)
    #x2 = np.array(xi2,dtype=ct.c_float)
    #x3 = np.array(xi3,dtype=ct.c_float)

    dsx = ct.c_float(Bsz/Nc)
    Ngb = ct.c_long(32)
    xc1 = ct.c_float(xic1)
    xc2 = ct.c_float(xic2)
    xc3 = ct.c_float(xic3)
    mass_particle = ct.c_float(pmass)
    posx1 = np.zeros((Nc,Nc),dtype=ct.c_float)
    posx2 = np.zeros((Nc,Nc),dtype=ct.c_float)
    sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    sps.cal_sph_sdens_arrays(x1,x2,x3,ct.c_float(Bsz),ct.c_long(Nc),dsx,Ngb,ct.c_long(Np),xc1,xc2,xc3,mass_particle,posx1,posx2,sdens)
    return sdens


sps.cal_sph_sdens_arrays_weights.argtypes =[np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     np.ctypeslib.ndpointer(dtype = ct.c_float),\
                                     ct.c_float,ct.c_long, \
                                     ct.c_float,ct.c_long,ct.c_long, \
                                     ct.c_float,ct.c_float,ct.c_float, \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float), \
                                     np.ctypeslib.ndpointer(dtype = ct.c_float)]
sps.cal_sph_sdens_arrays_weights.restype  = ct.c_int

def call_sph_sdens_arrays_weights(xi1,xi2,xi3,xic1,xic2,xic3,Bsz,Nc,Np,pmass):

    x1 = xi1.astype(ct.c_float)
    x2 = xi2.astype(ct.c_float)
    x3 = xi3.astype(ct.c_float)
    #x1 = np.array(xi1,dtype=ct.c_float)
    #x2 = np.array(xi2,dtype=ct.c_float)
    #x3 = np.array(xi3,dtype=ct.c_float)

    dsx = ct.c_float(Bsz/Nc)
    Ngb = ct.c_long(32)
    xc1 = ct.c_float(xic1)
    xc2 = ct.c_float(xic2)
    xc3 = ct.c_float(xic3)
    mass_particle = np.array(pmass,dtype=ct.c_float)
    posx1 = np.zeros((Nc,Nc),dtype=ct.c_float)
    posx2 = np.zeros((Nc,Nc),dtype=ct.c_float)
    sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    sps.cal_sph_sdens_arrays_weights(x1,x2,x3,ct.c_float(Bsz),ct.c_long(Nc),dsx,Ngb,ct.c_long(Np),xc1,xc2,xc3,mass_particle,posx1,posx2,sdens)
    return sdens
#--------------------------------------------------------------------
rsl = ct.CDLL("../improving_OM10_using_Illustris/lib_so_ray_tracing_single/lib/libraysingle.so")
rsl.ray_tracing_single.argtypes =[np.ctypeslib.ndpointer(dtype = ct.c_double),\
                                     np.ctypeslib.ndpointer(dtype = ct.c_double),\
                                     ct.c_int,ct.c_double,ct.c_double, \
                                     np.ctypeslib.ndpointer(dtype = ct.c_double)]
rsl.ray_tracing_single.restype  = ct.c_int

def call_ray_tracing_single(ai1,ai2,Nc,Bsz,zl):

    alpha1 = np.array(ai1,dtype=ct.c_double)
    alpha2 = np.array(ai2,dtype=ct.c_double)

    lensed_imgs_i = np.zeros((Nc,Nc),dtype=ct.c_double)

    rsl.ray_tracing_single(alpha1,alpha2,ct.c_int(Nc),ct.c_double(Bsz),ct.c_double(zl),lensed_imgs_i)
    return lensed_imgs_i
##---------------------------------------------------------------------------------
#scic = ct.CDLL("./libso_omp_cic/libcic.so")
#scic.cal_cic_sdens.argtypes =[np.ctypeslib.ndpointer(dtype = ct.c_float), \
                              #np.ctypeslib.ndpointer(dtype = ct.c_float), \
                              #ct.c_int,ct.c_float,ct.c_float,ct.c_float, \
                              #ct.c_int,ct.c_int, \
                              #np.ctypeslib.ndpointer(dtype = ct.c_float)]

#scic.cal_cic_sdens.restype  = ct.c_int

#def call_cic_sdens(x1,x2,Bsz,Nc,Np):

    #dsx = ct.c_float(Bsz/Nc)
    #xc1 = ct.c_float(0.0)
    #xc2 = ct.c_float(0.0)
    #x1 = np.array(x1,dtype=ct.c_float)
    #x2 = np.array(x2,dtype=ct.c_float)

    #sdens = np.zeros((Nc,Nc),dtype=ct.c_float)

    #scic.cal_cic_sdens(x1,x2,ct.c_int(Np),xc1,xc2,dsx,ct.c_int(Nc),ct.c_int(Nc),sdens)
    #return sdens
