#include <stdio.h>
#include <stdlib.h>
#include <math.h>

inline double deg2rad(double pha) {
	double res = 0;
	res = pha*M_PI/180.0;
	return res;
}

inline void lq_nie(double x1,double x2,double *lpar,double *alpha1,double *alpha2) {

    double xc1 = lpar[0];
    double xc2 = lpar[1];
    double q   = lpar[2];
    double rc  = lpar[3];
    double re  = lpar[4];
    double pha = lpar[5];

    double phirad = deg2rad(pha);
    double cosa = cos(phirad);
    double sina = sin(phirad);
	double phi,a1,a2;
	double xt1,xt2;

	xt1 = (x1-xc1)*cosa+(x2-xc2)*sina;
	xt2 = (x2-xc2)*cosa-(x1-xc1)*sina;
	phi = sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc);

	a1 = sqrt(q)/sqrt(1.0-q*q)*atan(sqrt(1.0-q*q)*xt1/(phi+rc/q));
	a2 = sqrt(q)/sqrt(1.0-q*q)*atanh(sqrt(1.0-q*q)*xt2/(phi+rc*q));

	*alpha1 = (a1*cosa-a2*sina)*re;
	*alpha2 = (a2*cosa+a1*sina)*re;
}

inline void tot_alphas(double x1, double x2,double *lpar, int npars, double *lpars, int nsubs, double *alpha1, double *alpha2) {

	int i,j;
    double al1, al2;
    lq_nie(x1,x2,lpar,&al1,&al2);

    double als1,als2;
    double * lpars_i = (double *)malloc(sizeof(double)*npars);
	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars;++j) {
			lpars_i[j] = lpars[i*npars+j];
		}
        lq_nie(x1,x2,lpars_i,&als1,&als2);
		al1 = al1+als1;
		al2 = al2+als2;
	}
	*alpha1 = al1;
	*alpha2 = al2;

    free(lpars_i);
}

inline void srcs_images(double xi1,double xi2,double *gpar,int npars,double *gpars,int nsubs,double *g_srcs) {
	int i,j;
	double g_srcs_tmp = 0.0;
    gauss_2d(xi1,xi2,gpar,&g_srcs_tmp);
    double * gpars_i = (double *)malloc(npars*sizeof(double));
    double g_lens_subs = 0;
	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars; ++j) {
			gpars_i[j] = gpars[i*npars+j];
		}
		gauss_2d(xi1,xi2,gpars_i,&g_lens_subs);
		g_srcs_tmp = g_srcs_tmp + g_lens_subs;
	}
	*g_srcs = g_srcs_tmp;
    free(gpars_i);
}


inline void single_ray_lensing(double xi1,double xi2,double * spar, int nspars, double * spars, int nssubs, double * lpar,int nlpars,double * lpars,int nlsubs,double *s_image,double *l_image){

	double al1,al2;
	double yi1,yi2;
	int i,k,l;

	tot_alphas(xi1,xi2,lpar,nlpars,lpars,nlsubs,&al1,&al2);

	yi1 = xi1-al1;
	yi2 = xi2-al2;

	//need to be improved
	double s_image_tmp;
	double l_image_tmp;
	srcs_images(xi1,xi2,spar,nspars,spars,nssubs,&s_image_tmp);
	srcs_images(yi1,yi2,spar,nspars,spars,nssubs,&l_image_tmp);
	*s_image = s_image_tmp;
	*l_image = l_image_tmp;
}

kernel void VectorAdd(
    global read_only double * xi1,
    global read_only double * xi2,
	global read_only double * spar,
	global read_only int nspars,
	global read_only double * spars,
	global read_only int nssubs,
	global read_only double * lpar,
	global read_only int nlpars,
	global read_only double * lpars,
	global read_only int nlsubs,
	global write_only double * s_image,
	global write_only double * l_image)
{
    int index = get_global_id(0);
	single_ray_lensing(xi1[index],xi2[index],spar,nspars,spars,nssubs,lpar,nlpars,lpars,nlsubs,&s_image[index],&l_image[index]);
}
