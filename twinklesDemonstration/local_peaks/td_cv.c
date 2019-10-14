#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double deg2rad(double pha) {
	double res = 0;
	res = pha*M_PI/180.0;
	return res;
}

void xy_rotate(double x1_in,double x2_in,double xc1,double xc2,double pha,double *x1_out,double *x2_out) {
    double phirad = deg2rad(pha);
	*x1_out = (x1_in - xc1)*cos(phirad)+(x2_in-xc2)*sin(phirad);
   	*x2_out = (x2_in - xc2)*cos(phirad)-(x1_in-xc1)*sin(phirad);
}

double gauss_2d(double x1,double x2,double *par) {
    //gpars = np.asarray([ylcs,xlcs,qls,aps,l_sigs,phis]);
	double x1new,x2new;
    xy_rotate(x1,x2,par[0],par[1],par[5],&x1new,&x2new);
	double res;
    res = par[3]*exp(-0.5*((x1new*x1new)*par[2]+(x2new*x2new)/par[2])/(par[4]*par[4]));
	return res;
}

double hfunc(double x1,double x2,double rcore, double qe) {
	double res;
    res = sqrt(qe*qe*(rcore*rcore+x1*x1)+x2*x2);
    return res;
}


double nie_phi(double x1,double x2,double re0,double rcore,double qe) {

	double res0;
    res0 = re0/sqrt(1-qe*qe);

	double al1,al2;
    al1 = res0*atan(x1*sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore));
    al2 = res0*atanh(x2*sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore*qe*qe));

	double res1;
    res1 = x1*al1+x2*al2;

	double res2;
    res2 = re0*rcore*log(sqrt(pow((hfunc(x1,x2,rcore,qe)+rcore),2.0)+(1-qe*qe)*x1*x1));

	double res;
    res = res1-res2;
    return res;
}

double nie_alphas(double x1,double x2,double re0,double rcore,double qe) {
	double res0;
    res0 = re0/sqrt(1-qe*qe);

	double al1,al2;
    al1 = atan(x1*sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore));
    al2 = atanh(x2*sqrt(1-qe*qe)/(hfunc(x1,x2,rcore,qe)+rcore*qe*qe));
    return res0*al1*res0*al2;
}

double nie_kappa(double x1,double x2,double re0,double rcore,double qe) {
	double res;
    res = re0/(2.0*sqrt(qe*qe*(rcore+x1*x1)+x2*x2));
    return res;
}

double nie_mu(double x1,double x2,double re0,double rcore,double qe) {
	double res;
    res = 1.0/(1.0-re0/hfunc(x1,x2,rcore,qe)-re0*re0*rcore/(hfunc(x1,x2,rcore,qe)*(pow(hfunc(x1,x2,rcore,qe)+rcore,2)+(1-qe*qe)*x1*x1)));
    return res;
}

double nie_td(double x1,double x2,double re0,double rcore,double qe,double ys1,double ys2) {

	double wx;
    wx = sqrt(qe*qe*(x1*x1+rcore*rcore)+x2*x2);

	double al1,al2;
    al1 = re0/sqrt(1-qe*qe)*atan(x1*sqrt(1-qe*qe)/(wx+rcore));
    al2 = re0/sqrt(1-qe*qe)*atanh(x2*sqrt(1-qe*qe)/(wx+qe*qe*rcore));

	double hx;
    hx = sqrt(pow((wx+rcore),2.0)+(1-qe*qe)*x1*x1);
	double phi;
    phi = x1*al1+x2*al2-re0*rcore*log(hx)+re0*qe*rcore*log((1+qe)*rcore);

    double Kc = 1.0;
    //Kc = (1.0+zl)/c*(Dl*Ds/Dls)
	double td;
    td = Kc*(0.5*(pow((x1-ys1),2.0)+pow((x2-ys2),2.0))-phi);

	return td;
}

void nie_all(double xi1,double xi2,double xc1,double xc2,double b,double s,double q,double rot,double ys1,double ys2,double *phi,double *td,double *al1,double *al2,double *kappa,double *mu,double *y1,double *y2) {

	double x1,x2;
    xy_rotate(xi1,xi2,xc1,xc2,rot,&x1,&x2);

	double wx;
    wx = sqrt(q*q*(x1*x1+s*s)+x2*x2);

    *al1 = b/sqrt(1-q*q)*atan(x1*sqrt(1-q*q)/(wx+s));
    *al2 = b/sqrt(1-q*q)*atanh(x2*sqrt(1-q*q)/(wx+q*q*s));

    *kappa = b/(2.0*wx);

	double hx;
    hx = sqrt(pow((wx+s),2.0)+(1-q*q)*x1*x1);
    *phi = x1*(*al1)+x2*(*al2)-b*s*log(hx)+b*q*s*log((1+q)*s);

    double Kc = 1.0;
    //Kc = (1.0+zl)/c*(Dl*Ds/Dls)
    *td = Kc*(0.5*(pow((x1-ys1),2.0)+pow((x2-ys2),2.0))-(*phi));

	double yt1,yt2;
    yt1 = x1-(*al1);
    yt2 = x2-(*al2);

    xy_rotate(yt1,yt2,xc1,xc2,-rot,y1,y2);

	double demon1,demon2;
    demon1 = (pow((wx+s),2)+(1.0-q*q)*x1*x1)*wx;
    demon2 = ((pow((wx+q*q*s),2)-(1.0-q*q)*x2*x2)*wx);
	double y11,y22,y12,y21;
    y11 = 1-b*(wx*(wx+s)-q*q*x1*x1)/demon1;
    y22 = 1-b*(wx*(wx+q*q*s)-x2*x2)/demon2;
    y12 = -b*x1*x2/demon1;
    y21 = -b*x1*x2*q*q/demon2;

    *mu = 1.0/(y11*y22-y12*y21);
}

void arbitrary_td(double *x1,double *x2,double *phi,int nx,int ny,double ys1,double ys2,double *td) {

    double Kc = 1.0;
    //Kc = (1.0+zl)/c*(Dl*Ds/Dls)
	int i,j,index;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			index = i*ny+j;
			td[index] = Kc*(0.5*(pow((x1[index]-ys1),2.0)+pow((x2[index]-ys2),2.0))-phi[index]);
		}
	}
}

void arbitrary_td_map(double *al1,double *al2,double *phi,int nx,int ny,double *td) {

    double Kc = 1.0;
    //Kc = (1.0+zl)/c*(Dl*Ds/Dls)
	int i,j,index;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			index = i*ny+j;
			td[index] = Kc*(0.5*(pow(al1[index],2.0)+pow(al2[index],2.0))-phi[index]);
		}
	}
}
