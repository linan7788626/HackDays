#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <omp.h>

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

void forward_cic(float *cic_in,float *x_in,float *y_in,float bsx,float bsy,int nx,int ny,int np,float *cic_out) {
    float dx = bsx/nx;
    float dy = bsy/ny;
    float xc = bsx/2.0;
    float yc = bsy/2.0;
    float wx,wy;
    float xp,yp,zp;

    int i;
    int ip,jp;

    for (i=0;i<np;i++) {
        xp = (x_in[i]+xc)/dx-0.5;
        yp = (y_in[i]+yc)/dy-0.5;
		zp = cic_in[i];

        ip = (int)xp;
        jp = (int)yp;

		if (ip<0||ip>(nx-2)||jp<0||jp>(ny-2)) continue;
        wx = 1.0-(xp-(float)ip);
        wy = 1.0-(yp-(float)jp);

        cic_out[ip*ny+jp] += wx*wy*zp;
        cic_out[ip*ny+(jp+1)] += wx*(1.0-wy)*zp;
        cic_out[(ip+1)*ny+jp] += (1.0-wx)*wy*zp;
        cic_out[(ip+1)*ny+(jp+1)] += (1.0-wx)*(1.0-wy)*zp;
    }
}

//--------------------------------------------------------------------
void lanczos_diff_2_tag(float *m1, float *m2, float *m11, float *m12, float *m21, float *m22, float Dcell, int Ncc, int dif_tag) {
    int i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
    int index;

    for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {

		//i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;
		//i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;
		//j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;
		//j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;

		//if (j_m3<0||j_p3>(Ncc-1)||i_m3<0||i_p3>(Ncc-1)) continue;
		//if (j_m2<0||j_p2>(Ncc-1)||i_m2<0||i_p2>(Ncc-1)) continue;
		//if (j_m1<0||j_p1>(Ncc-1)||i_m1<0||i_p1>(Ncc-1)) continue;

        if (i==0) {i_m1 = Ncc-1;i_m2 = Ncc-2;i_m3 = Ncc-3;}
        else if (i==1) {i_m1 = 0;i_m2 = Ncc-1;i_m3 = Ncc-2;}
        else if (i==2) {i_m1 = 1;i_m2 = 0;i_m3 = Ncc-1;}
        else {i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;}
        if (j==0) {j_m1 = Ncc-1;j_m2 = Ncc-2;j_m3 = Ncc-3;}
        else if (j==1) {j_m1 = 0;j_m2 = Ncc-1;j_m3 = Ncc-2;}
        else if (j==2) {j_m1 = 1;j_m2 = 0;j_m3 = Ncc-1;}
        else {j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;}
        if (i==Ncc-1) {i_p1 = 0;i_p2 = 1;i_p3 = 2;}
        else if (i==Ncc-2) {i_p1 = Ncc-1;i_p2 = 0;i_p3 = 1;}
        else if (i==Ncc-3) {i_p1 = Ncc-2;i_p2 = Ncc-1;i_p3 = 0;}
        else {i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;}
        if (j==Ncc-1) {j_p1 = 0;j_p2 = 1;j_p3 = 2;}
        else if (j==Ncc-2) {j_p1 = Ncc-1;j_p2 = 0;j_p3 = 1;}
        else if (j==Ncc-2) {j_p1 = Ncc-2;j_p2 = Ncc-1;j_p3 = 0;}
        else {j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;}

        index = i*Ncc+j;
        if (dif_tag==-1) {
            m11[index] = (m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])/(2.0*Dcell);
            m12[index] = (m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])/(2.0*Dcell);
            m21[index] = (m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])/(2.0*Dcell);
            m22[index] = (m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])/(2.0*Dcell);
        }

        if (dif_tag==0) {
            m11[index] = (m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])*2.0/3.0/Dcell
            - (m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])/12.0/Dcell;
            m22[index] = (m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])*2.0/3.0/Dcell
            - (m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])/12.0/Dcell;
            m21[index] = (m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])*2.0/3.0/Dcell
            - (m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])/12.0/Dcell;
            m12[index] = (m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])*2.0/3.0/Dcell
            - (m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])/12.0/Dcell;
        }

        if (dif_tag==1) {
            m11[index] =(1.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
                         + 2.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
                         + 3.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(28.0*Dcell);
            m22[index] =(1.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
                         + 2.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
                         + 3.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(28.0*Dcell);
            m12[index] =(1.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
                         + 2.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
                         + 3.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(28.0*Dcell);
            m21[index] =(1.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
                         + 2.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
                         + 3.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(28.0*Dcell);
        }

        if (dif_tag==2) {
            m11[index] = (5.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
                          + 4.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
                          + 1.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(32.0*Dcell);
            m22[index] = (5.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
                          + 4.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
                          + 1.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(32.0*Dcell);
            m12[index] = (5.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
                          + 4.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
                          + 1.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(32.0*Dcell);
            m21[index] = (5.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
                          + 4.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
                          + 1.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(32.0*Dcell);
        }

        if (dif_tag==3) {
            m11[index] = (58.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
                          + 67.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
                          + 22.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(252.0*Dcell);
            m22[index] = (58.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
                          + 67.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
                          - 22.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(252.0*Dcell);
            m12[index] = (58.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
                          + 67.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
                          - 22.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(252.0*Dcell);
            m21[index] = (58.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
                          + 67.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
                          - 22.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(252.0*Dcell);
        }
    }
}

void xy_rotate(float x1_in,float x2_in,float xc1,float xc2,float pha,float *x1_out,float *x2_out) {
    float phirad = deg2rad(pha);
	*x1_out = (x1_in - xc1)*cos(phirad)+(x2_in-xc2)*sin(phirad);
	*x2_out = (x2_in - xc2)*cos(phirad)-(x1_in-xc1)*sin(phirad);
}

void gauss_2d(float x1,float x2,float *par,float *res) {
    //gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis]);
	float x1new,x2new;
    xy_rotate(x1,x2,par[0], par[1],par[5],&x1new,&x2new);
	float re_eff = (x1new*x1new)*par[2]+(x2new*x2new)/par[2];
	if (re_eff>4.0) {
		*res = 0.0;
	}
	else {
		*res = par[3]*exp(-0.5*(re_eff)/(par[4]*par[4]));

	}
}

void tophat_2d(float x1,float x2,float *par,float *res) {
	float x1new,x2new;
    xy_rotate(x1,x2,par[0],par[1], par[5],&x1new,&x2new);
	float r_ell;
	r_ell = sqrt((x1new*x1new)*par[2]+(x2new*x2new)/par[2]);
	if (r_ell>=par[4]) {
		*res = -1.0;
	}
	else {
		*res = 10000.0;
	}
}

void lq_nie(float x1,float x2,float *lpar,float *alpha1,float *alpha2) {

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
	float xt1,xt2;

	xt1 = (x1-xc1)*cosa+(x2-xc2)*sina;
	xt2 = (x2-xc2)*cosa-(x1-xc1)*sina;
	phi = sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc);

	a1 = sqrt(q)/sqrt(1.0-q*q)*atan(sqrt(1.0-q*q)*xt1/(phi+rc/q));
	a2 = sqrt(q)/sqrt(1.0-q*q)*atanh(sqrt(1.0-q*q)*xt2/(phi+rc*q));

	*alpha1 = (a1*cosa-a2*sina)*re;
	*alpha2 = (a2*cosa+a1*sina)*re;
}

void find_critical_curve(float *mu,int nx,int ny,float* res) {

	int i,j,index,sign_t=0;
	int im1,ip1,jm1,jp1;
	for (i = 0; i < nx; ++i) for (j = 0; j < ny; ++j) {
		index = i*ny+j;
		//im1 = i-1;
		//ip1 = i+1;
		//jm1 = j-1;
		//jp1 = j+1;

		//if (im1<0||jm1<0||ip1>(nx-1)||jp1>(ny-1)) continue;

        if (i==0) {im1 = nx-1;}
        else {im1 = i-1;}
        if (j==0) {jm1 = ny-1;}
        else {jm1 = j-1;}
        if (i==nx-1) {ip1 = 0;}
        else {ip1 = i+1;}
        if (j==ny-1) {jp1 = 0;}
        else {jp1 = j+1;}

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
void tot_lq(float x1, float x2,float *lpar, int npars, float *lpars, int nsubs, float *y1, float *y2) {

	int i,j;
    float al1, al2;
    lq_nie(x1,x2,lpar,&al1,&al2);

    float als1,als2;
    float * lpars_i = (float *)malloc(sizeof(float)*npars);
	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars;++j) {
			lpars_i[j] = lpars[i*npars+j];
		}
        lq_nie(x1,x2,lpars_i,&als1,&als2);
		al1 = al1+als1;
		al2 = al2+als2;
	}

	*y1 = x1-al1;
	*y2 = x2-al2;
    free(lpars_i);
}
void tot_alphas(float x1, float x2,float *lpar, int npars, float *lpars, int nsubs, float *alpha1, float *alpha2) {

	int i,j;
    float al1, al2;
    lq_nie(x1,x2,lpar,&al1,&al2);

    float als1,als2;
    float * lpars_i = (float *)malloc(sizeof(float)*npars);
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

void refine_critical(float * xi1,float * xi2,int nx1,int nx2,float * lpar,int npars,float * lpars, int nsubs,float * critical,int clen, int nfiner, float * yi1,float *yi2) {

	int i,j,k=0,m,n,index;
    float dsx = xi1[nx2+1]-xi1[0];
    float dsf = dsx/nfiner;

    float * xt1 = (float *)malloc(clen*nfiner*nfiner*sizeof(float));
    float * xt2 = (float *)malloc(clen*nfiner*nfiner*sizeof(float));

	for (i = 0; i < nx1; ++i) for (j = 0; j < nx2; ++j){
		index = i*nx2+j;
		if (critical[index]>0) {
			for (m = 0; m < nfiner; ++m) for (n = 0; n < nfiner; ++n){
        	    xt1[k*nfiner*nfiner+m*nfiner+n] = xi1[index]+(dsf*(1-nfiner)*0.5)+dsf*m;
        	    xt2[k*nfiner*nfiner+m*nfiner+n] = xi2[index]+(dsf*(1-nfiner)*0.5)+dsf*n;
			}
			k = k+1;
		}
	}
	for (i = 0; i < clen*nfiner*nfiner; i++) {
		tot_lq(xt1[i],xt2[i],lpar,npars,lpars,nsubs,&yi1[i],&yi2[i]);
	}
	free(xt1);
	free(xt2);
}
void srcs_images(float xi1,float xi2,float *gpar,int npars,float *gpars,int nsubs,float *g_srcs) {
	int i,j;
	float g_srcs_tmp = 0.0;
    gauss_2d(xi1,xi2,gpar,&g_srcs_tmp);
    float * gpars_i = (float *)malloc(npars*sizeof(float));
    float g_lens_subs = 0;
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

void lens_images(float *xi1,float *xi2,int nx1,int nx2,float *gpar,int npars,float *gpars,int nsubs,float *g_lens) {
	int i,j,k,l,index;
	for (k = 0; k < nx1; ++k) for (l = 0; l < nx2; ++l){
		index = k*nx2+l;
		gauss_2d(xi1[index],xi2[index],gpar,&g_lens[index]);
	}
    float * gpars_i = (float *)malloc(npars*sizeof(float));
    float * g_lens_subs = (float *)malloc(nx1*nx2*sizeof(float));
	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars; ++j) {
			gpars_i[j] = gpars[i*npars+j];
		}
		for (k = 0; k < nx1; ++k) for (l = 0; l < nx2; ++l){
			index = k*nx2+l;
			gauss_2d(xi1[i],xi2[i],gpars_i,&g_lens_subs[index]);
			g_lens[index] = g_lens[index] + g_lens_subs[index];
		}
	}
    free(gpars_i);
    free(g_lens_subs);
}
void mmbr_images(float *xi1,float *xi2,int nx1,int nx2,float *gpar,int npars,float *gpars,int nsubs,float *g_edge) {

	int i,j,k,l,index;
    float * g_lens = (float *)malloc(nx1*nx2*sizeof(float));
	for (i = 0; i < nx1*nx2; i++) {
		tophat_2d(xi1[i],xi2[i],gpar,&g_lens[i]);
	}
    find_critical_curve(g_lens,nx1,nx2,g_edge);

    float * gpars_i = (float *)malloc(npars*sizeof(float));
    float * g_lens_subs = (float *)malloc(nx1*nx2*sizeof(float));
    float * g_edge_subs = (float *)malloc(nx1*nx2*sizeof(float));

	for (i = 0; i < nsubs; ++i) {
		for (j = 0; j < npars; ++j) {
			gpars_i[j] = gpars[i*npars+j];
		}
		for (i = 0; i < nx1*nx2; i++) {
			tophat_2d(xi1[i],xi2[i],gpars_i,&g_lens_subs[i]);
		}
        find_critical_curve(g_lens_subs,nx1,nx2,g_edge_subs);
		for (k = 0; k < nx1; ++k) for (l = 0; l < nx2; ++l){
			index = k*nx2+l;
			g_edge[index] = g_edge[index] + g_edge_subs[index];
		}
	}
    free(gpars_i);
    free(g_lens_subs);

	for (k = 0; k < nx1; ++k) for (l = 0; l < nx2; ++l){
		index = k*nx2+l;
		if (g_edge[index]>0.0) {
			g_edge[index] = 1.0;
		}
	}
}

void find_caustics(float *xi1,float *xi2,int nx1,int nx2,float dsx,float *critical,float *lpar,int nlpars,float *lpars,int nlsubs,float *caustic) {

	int i,j,index;
    int clen = 0;
	int nfiner = 4;

	for (i = 0; i < nx1; ++i) for (j = 0; j < nx2; ++j){
		index = i*nx2+j;
		if (critical[index]>0) {
			clen = clen+1;
		}
	}
	int ylen = clen*nfiner*nfiner;

    float * yif1 = (float *)malloc(ylen*sizeof(float));
    float * yif2 = (float *)malloc(ylen*sizeof(float));

	refine_critical(xi1,xi2,nx1,nx2,lpar,nlpars,lpars,nlsubs,critical,clen,nfiner,yif1,yif2);

    float bsz;
    bsz = dsx*nx1;
    float * img_in = (float *)malloc(ylen*sizeof(float));

	for (i = 0; i < ylen; ++i) {
		img_in[i] = 1.0;
	}
    forward_cic(img_in,yif1,yif2,bsz,bsz,nx1,nx2,ylen,caustic);

	free(yif1);
	free(yif2);
	free(img_in);

	for (i = 0; i < nx1; ++i) for (j = 0; j < nx2; ++j){
		index = i*nx2+j;
		if (caustic[index]>0) {
			caustic[index] = 1;
		}
	}
}


void all_about_lensing(float *xi1,float *xi2,int nx1,int nx2,float * spar, int nspars, float * spars, int nssubs, float * lpar,int nlpars,float * lpars,int nlsubs,float *s_image,float *g_lensimage,float *critical,float *caustic){

    float * al1 = (float *)malloc(sizeof(float)*nx1*nx2);
    float * al2 = (float *)malloc(sizeof(float)*nx1*nx2);
	float yi1,yi2;
	int i,k,l;

#pragma omp parallel num_threads(4)    \
    shared(xi1,xi2,nx1,nx2,lpar,nlpars,lpars,nlsubs,spar,nspars,spars,nssubs,s_image,g_lensimage) \
    private(i,yi1,yi2)
    {
    #pragma omp for schedule(dynamic,16)
		for (i = 0; i < nx1*nx2; i++) {
			tot_alphas(xi1[i],xi2[i],lpar,nlpars,lpars,nlsubs,&al1[i],&al2[i]);

			yi1 = xi1[i]-al1[i];
			yi2 = xi2[i]-al2[i];
			//gauss_2d(xi1[i],xi2[i],spar,&s_image[i]);
			//srcs_images(xi1[i],xi2[i],spar,nspars,spars,nssubs,&s_image[i]);
			//srcs_images(yi1,yi2,spar,nspars,spars,nssubs,&g_lensimage[i]);
		}
    }
//------------------------------------------------------------------------
    float dsx = xi1[nx2+1]-xi1[0];
    float * a11 = (float *)malloc(nx1*nx2*sizeof(float));
    float * a12 = (float *)malloc(nx1*nx2*sizeof(float));
    float * a21 = (float *)malloc(nx1*nx2*sizeof(float));
    float * a22 = (float *)malloc(nx1*nx2*sizeof(float));

	lanczos_diff_2_tag(al1,al2,a21,a22,a11,a12,dsx,nx1,-1);
    free(al1);
    free(al2);

    float * imu = (float *)malloc(nx1*nx2*sizeof(float));
	int index;
	for (k = 0; k < nx1; k++) for (l = 0; l < nx2; l++) {
		index = k*nx2+l;
		imu[index] = (1.0-(a11[index]+a22[index])+a11[index]*a22[index]-a12[index]*a21[index]);
	}

	free(a11);
	free(a12);
	free(a21);
	free(a22);

    find_critical_curve(imu,nx1,nx2,critical);
	free(imu);
	find_caustics(xi1,xi2,nx1,nx2,dsx,critical,lpar,nlpars,lpars,nlsubs,caustic);
}

void single_ray_lensing(float xi1,float xi2,float * spar, int nspars, float * spars, int nssubs, float * lpar,int nlpars,float * lpars,int nlsubs,float *s_image,float *l_image){

	float al1,al2;
	float yi1,yi2;
	int i,k,l;

	tot_alphas(xi1,xi2,lpar,nlpars,lpars,nlsubs,&al1,&al2);

	yi1 = xi1-al1;
	yi2 = xi2-al2;

	//need to be improved
	float s_image_tmp;
	float l_image_tmp;
	srcs_images(xi1,xi2,spar,nspars,spars,nssubs,&s_image_tmp);
	srcs_images(yi1,yi2,spar,nspars,spars,nssubs,&l_image_tmp);
	*s_image = s_image_tmp;
	*l_image = l_image_tmp;
}
//------------------------------------------------------------------------
void cal_cc(float *xi1,float *xi2,float *al1,float *al2,int nx1,int nx2,float *lpar,int nlpars,float *lpars,int nlsubs,float *critical,float *caustic){
	int i,j,k,l;
    float dsx = xi1[nx2+1]-xi1[0];
    float * a11 = (float *)malloc(nx1*nx2*sizeof(float));
    float * a12 = (float *)malloc(nx1*nx2*sizeof(float));
    float * a21 = (float *)malloc(nx1*nx2*sizeof(float));
    float * a22 = (float *)malloc(nx1*nx2*sizeof(float));

	lanczos_diff_2_tag(al1,al2,a21,a22,a11,a12,dsx,nx1,-1);
    free(al1);
    free(al2);

    float * imu = (float *)malloc(nx1*nx2*sizeof(float));
	int index;
	for (k = 0; k < nx1; k++) for (l = 0; l < nx2; l++) {
		index = k*nx2+l;
		imu[index] = (1.0-(a11[index]+a22[index])+a11[index]*a22[index]-a12[index]*a21[index]);
	}

	free(a11);
	free(a12);
	free(a21);
	free(a22);

    find_critical_curve(imu,nx1,nx2,critical);
	free(imu);
	find_caustics(xi1,xi2,nx1,nx2,dsx,critical,lpar,nlpars,lpars,nlsubs,caustic);
}
