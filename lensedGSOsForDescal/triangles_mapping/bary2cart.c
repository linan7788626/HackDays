#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef enum { false, true } bool;

double signP(double *p0, double *p1, double *p2) {
	return (p0[0] - p2[0]) * (p1[1] - p2[1]) - (p1[0] - p2[0]) * (p0[1] - p2[1]);
}

bool PIT(double *pt, double *v0, double *v1, double *v2) {
	bool b0, b1, b2;

	b0 = signP(pt, v0, v1) < 0.0f;
	b1 = signP(pt, v1, v2) < 0.0f;
	b2 = signP(pt, v2, v0) < 0.0f;

	return ((b0 == b1) && (b1 == b2));
}

double triAera(double * p0, double *p1, double *p2) {
	double res = 0.0;

	/*res = a[0][0]*((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) -a[0][1]*(a[1][0]*a[2][2] - a[2][0]*a[1][2]) + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1])*/
	res = p0[0]*(p1[1]-p2[1])-p0[1]*(p1[0]-p2[0])+(p1[0]*p2[1]-p2[0]*p1[1]);
	return res;
}

void bary2cart(double *p0, double *p1, double *p2, double *bary, double *pt) {
    pt[0] = bary[0] * p0[0] + bary[1] * p1[0] + bary[2] * p2[0];
    pt[1] = bary[0] * p0[1] + bary[1] * p1[1] + bary[2] * p2[1];
}

double dotP(double *v, double *u, int n){
	int i;
	double res = 0.0;
	for (i = 0; i < n; i++) {
		res += v[i]*u[i];
	}
	return res;
}

void Cart2Bary(double *p, double *a, double *b, double *c, double *bary) {
    /*double *v0 = b - a;*/
	/*double *v1 = c - a;*/
	/*double *v2 = p - a;*/

	double *v0 = (double *)malloc(sizeof(double)*2);
	v0[0] = b[0]-a[0];
	v0[1] = b[1]-a[1];
	double *v1 = (double *)malloc(sizeof(double)*2);
	v1[0] = c[0]-a[0];
	v1[1] = c[1]-a[1];
	double *v2 = (double *)malloc(sizeof(double)*2);
	v2[0] = p[0]-a[0];
	v2[1] = p[1]-a[1];

    double d00 = dotP(v0, v0, 2);
    double d01 = dotP(v0, v1, 2);
    double d11 = dotP(v1, v1, 2);
    double d20 = dotP(v2, v0, 2);
    double d21 = dotP(v2, v1, 2);
    double denom = d00 * d11 - d01 * d01;

    bary[1] = (d11 * d20 - d01 * d21) / denom;
    bary[2] = (d00 * d21 - d01 * d20) / denom;
    bary[0] = 1.0f - bary[1] - bary[2];
}

// p0, p1 and p2  and the points that make up this triangle
void cart2bary(double *pt, double *p0, double *p1, double *p2, double *bary) {
    double tri_area = triAera(p0,p1,p2);
    bary[0] = triAera(pt,p1,p2) / tri_area;
    bary[1] = triAera(pt,p2,p0) / tri_area;
    bary[2] = triAera(pt,p0,p1) / tri_area;
}

//--------------------------------------------------------------------
void tria_vectors(int i,int j,double *xgrids1,double *xgrids2,int nc,double *v0,double *v1,double *v2) {
    int ip1 = i+1;
    int jp1 = j+1;

    v0[0] = xgrids1[i*nc+j];
	v0[1] = xgrids2[i*nc+j];
    v1[0] = xgrids1[ip1*nc+j];
	v1[1] = xgrids2[ip1*nc+j];
    v2[0] = xgrids1[ip1*nc+jp1];
	v2[1] = xgrids2[ip1*nc+jp1];
}



void trib_vectors(int i,int j,double *xgrids1,double *xgrids2,int nc,double *v0,double *v1,double *v2) {

    int ip1 = i+1;
    int jp1 = j+1;

    v0[0] = xgrids1[i*nc+j];
	v0[1] = xgrids2[i*nc+j];
    v1[0] = xgrids1[ip1*nc+jp1];
	v1[1] = xgrids2[ip1*nc+jp1];
    v2[0] = xgrids1[i*nc+jp1];
	v2[1] = xgrids2[i*nc+jp1];
}


void mapping_triangles(double *ys,double *xgrids1,double *xgrids2,double *lgrids1,double *lgrids2, int nc, double *xroots) {

	int i,j;
    int ntris = nc-1;
    int ncount = 0;

	double *xv0 = (double *)malloc(sizeof(double)*2);
	double *xv1 = (double *)malloc(sizeof(double)*2);
	double *xv2 = (double *)malloc(sizeof(double)*2);

	double *yv0 = (double *)malloc(sizeof(double)*2);
	double *yv1 = (double *)malloc(sizeof(double)*2);
	double *yv2 = (double *)malloc(sizeof(double)*2);

	double *baryc = (double *)malloc(sizeof(double)*3);

	double *xroots_tmp = (double *)malloc(sizeof(double)*2);

	for (i = 0; i < ntris; i++) {
		for (j = 0; j < ntris; j++) {
            tria_vectors(i,j,xgrids1,xgrids2,nc,xv0,xv1,xv2);
            tria_vectors(i,j,lgrids1,lgrids2,nc,yv0,yv1,yv2);
            if (PIT(ys,yv0,yv1,yv2)){
                cart2bary(ys,yv0,yv1,yv2,baryc);
                bary2cart(xv0,xv1,xv2,baryc,xroots_tmp);
				xroots[ncount*2+0] = xroots_tmp[0];
				xroots[ncount*2+1] = xroots_tmp[1];
                ncount = ncount + 1;
			}
		}
	}

	for (i = 0; i < ntris; i++) {
		for (j = 0; j < ntris; j++) {
            trib_vectors(i,j,xgrids1,xgrids2,nc,xv0,xv1,xv2);
            trib_vectors(i,j,lgrids1,lgrids2,nc,yv0,yv1,yv2);
            if (PIT(ys,yv0,yv1,yv2)){
                cart2bary(ys,yv0,yv1,yv2,baryc);
                bary2cart(xv0,xv1,xv2,baryc,xroots_tmp);
				xroots[ncount*2+0] = xroots_tmp[0];
				xroots[ncount*2+1] = xroots_tmp[1];
                ncount = ncount + 1;
			}
		}
	}

	free(xv0);
	free(xv1);
	free(xv2);
	free(yv0);
	free(yv1);
	free(yv2);
	free(baryc);
	free(xroots_tmp);
}
//--------------------------------------------------------------------

int main() {
	double * pt = (double *)malloc(sizeof(double)*2);
	double * v0 = (double *)malloc(sizeof(double)*2);
	double * v1 = (double *)malloc(sizeof(double)*2);
	double * v2 = (double *)malloc(sizeof(double)*2);

	pt[0] = 0.3;
	pt[1] = 0.3;

	v0[0] = 0.0;
	v0[1] = 0.0;

	v1[0] = 1.0;
	v1[1] = 0.0;

	v2[0] = 0.0;
	v2[1] = 1.0;

	bool res;

	res = PIT(pt,v0,v1,v2);
	printf("%s \n", (res)?"True":"False");

	if (res) {
		double * barv = (double *)malloc(sizeof(double)*3);
		/*Cart2Bary(pt, v0, v1, v2, barv);*/
		cart2bary(pt, v0, v1, v2, barv);
		printf("%lf %lf %lf\n", barv[0], barv[1], barv[2]);

		double * v02 = (double *)malloc(sizeof(double)*2);
		double * v12 = (double *)malloc(sizeof(double)*2);
		double * v22 = (double *)malloc(sizeof(double)*2);

		v02[0] = 0.0*3.0;
		v02[1] = 0.0*3.0;

		v12[0] = 1.0*3.0;
		v12[1] = 0.0*3.0;

		v22[0] = 0.0*3.0;
		v22[1] = 1.0*3.0;

		double * pt2 = (double *)malloc(sizeof(double)*2);
		bary2cart(v02, v12, v22, barv, pt2);
		printf("%lf %lf\n", pt2[0], pt2[1]);

		free(v02);
		free(v12);
		free(v22);
		free(barv);
		free(pt2);
	}

	free(pt);
	free(v0);
	free(v1);
	free(v2);

	return 0;
}
