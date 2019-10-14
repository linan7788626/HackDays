#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int sign(float x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

float deg2rad(float pha) {
	float res = 0;
	res = pha*M_PI/180.0;
	return res;
}

void xy_rotate(float *x1_in,float *x2_in,int nx1,int nx2,float xc1,float xc2,float pha,float *x1_out,float *x2_out) {
    float phirad = deg2rad(pha);
	int i,j,index;
	for (i=0;i<nx1;++i) for (j=0;j<nx2;++j){
		index = i*nx2+j;
		x1_out[index] = (x1_in[index] - xc1)*cos(phirad)+(x2_in[index]-xc2)*sin(phirad);
    	x2_out[index] = (x2_in[index] - xc2)*cos(phirad)-(x1_in[index]-xc1)*sin(phirad);
	}
}

void gauss_2d(float *x1,float *x2,int nx1,int nx2,float *par,float *res) {
    //gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis]);
	float *x1new = (float *)malloc(nx1*nx2*sizeof(float));
	float *x2new = (float *)malloc(nx1*nx2*sizeof(float));
    xy_rotate(x1,x2,nx1,nx2,par[0], par[1],par[5],x1new,x2new);
	int i,j,index;
	for (i=0;i<nx1;++i) for (j=0;j<nx2;++j) {
		index = i*nx2+j;
    	res[index] = par[3]*exp(-0.5*((x1new[index]*x1new[index])*par[2]+(x2new[index]*x2new[index])/par[2])/(par[4]*par[4]));
	}
	free(x1new);
	free(x2new);
}

void tophat_2d(float *x1,float *x2, int nx1, int nx2,float *par,float *res) {
    //gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis]);
	float *x1new = (float *)malloc(nx1*nx2*sizeof(float));
	float *x2new = (float *)malloc(nx1*nx2*sizeof(float));
    xy_rotate(x1,x2,nx1,nx2,par[0], par[1], par[5],x1new,x2new);
	int i,j,index;
	float r_ell;
	for (i=0;i<nx1;++i) for (j=0;j<nx2;++j){
		index = i*nx2+j;
    	r_ell = sqrt((x1new[index]*x1new[index])*par[2]+(x2new[index]*x2new[index])/par[2]);
		if (r_ell>=par[4]) {
			res[index] = -1.0;
		}
		else {
			res[index] = 10000.0;
		}
	}
	free(x1new);
	free(x2new);
}

void lq_nie(float *x1,float *x2,int nx1,int nx2,float *lpar,float *alpha1,float *alpha2) {

    float xc1 = lpar[0];
    float xc2 = lpar[1];
    float q   = lpar[2];
    float rc  = lpar[3];
    float re  = lpar[4];
    float pha = lpar[5];

    float phirad = deg2rad(pha);
    float cosa = cos(phirad);
    float sina = sin(phirad);
	float phi,a1,a2;
	float *xt1 = (float *)malloc(sizeof(float)*nx1*nx2);
	float *xt2 = (float *)malloc(sizeof(float)*nx1*nx2);

	int i,j,index;
	for (i = 0;i<nx1;++i) for (j = 0;j<nx2;++j) {
		index = i*nx2+j;
		xt1[index] = (x1[index]-xc1)*cosa+(x2[index]-xc2)*sina;
    	xt2[index] = (x2[index]-xc2)*cosa-(x1[index]-xc1)*sina;
		phi = sqrt(xt2[index]*xt2[index]+xt1[index]*q*xt1[index]*q+rc*rc);

    	a1 = sqrt(q)/sqrt(1.0-q*q)*atan(sqrt(1.0-q*q)*xt1[index]/(phi+rc/q));
    	a2 = sqrt(q)/sqrt(1.0-q*q)*atanh(sqrt(1.0-q*q)*xt2[index]/(phi+rc*q));

    	alpha1[index] = (a1*cosa-a2*sina)*re;
    	alpha2[index] = (a2*cosa+a1*sina)*re;
	}
	free(xt1);
	free(xt2);
}

//void refine_critical(float *critical,float *xi1,float *xi2,float *ai1,float *ai2,int n,float dsx,int nfiner,float *yi1,float *yi2) {
//	int i,j;
//	int ncount = 0;
//	for (i = 0; i < n; ++i) {
//		if (critical[i]>0) {
//			ncount = ncount+1;
//		}
//	}
//	float *x1tmp = (float *) malloc(sizeof(float)*ncount);
//	float *x2tmp = (float *) malloc(sizeof(float)*ncount);
//	ncount = 0;
//	for (i = 0; i < n; ++i) {
//		if (critical[i]>0) {
//			x1tmp[ncount] = xi1[i];
//			x2tmp[ncount] = xi2[i];
//			ncount = ncount+1;
//		}
//	}
//
//	int k;
//    int dsf = dsx/nfiner/2;
//	float x1t,x2t;
//	for (k = 0; k<ncount; ++k) {
//		x1t = x1tmp[k];
//		x2t = x2tmp[k];
//		for (i = 0; i<nfiner; ++i) for (j = 0; j<nfiner; ++j){
//			x1t = x1tmp[k]+dsf*(1.0-nfiner)*0.5+dsf*i;
//			x2t = x2tmp[k]+dsf*(1.0-nfiner)*0.5+dsf*i;
//		}
//	}
//    for i in xrange(nfiner):
//        for j in xrange(nfiner):
//            x1tmp = xi1[critical>0]+(dsf*(1-nfiner)*0.5)+dsf*i
//            x2tmp = xi2[critical>0]+(dsf*(1-nfiner)*0.5)+dsf*j
//
//            yift1[:,i,j],yift2[:,i,j] = source_plane_finer(x1tmp,x2tmp,lpar,lpars)
//
//    return yift1,yift2
//}

void find_critical_curve(float *mu,int nx,int ny,float* res) {

	int i,j,index,sign_t=0;
	int im1,ip1,jm1,jp1;
	for (i = 0; i < nx; ++i) for (j = 0; j < ny; ++j) {
		index = i*ny+j;
		im1 = i-1;
		ip1 = i+1;
		jm1 = j-1;
		jp1 = j+1;

		if (im1<0||jm1<0||ip1>(nx-1)||jp1>(ny-1)) continue;

		sign_t = sign(mu[index])*(sign(mu[im1*ny+j])
								 +sign(mu[i*ny+jm1])
								 +sign(mu[ip1*ny+j])
								 +sign(mu[i*ny+jp1]));
		if (sign_t < 4) {
			res[index] = 1.0;
		}
		else {
			res[index] = 0.0;
		}
	}
}

//int main(int argc, const char *argv[])
//{
//
//	return 0;
//}

