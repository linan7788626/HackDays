/*
 * find the critical curves and caustics for an isothermal elliptical density distribution plus point lenses
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#define  MIN(a,b)       (((a) <= (b)) ? (a) : (b))
#define  TINY	1.0e-20

int nlenses=2;			/* number of point lenses */
double *xl, *yl, *ml;		/* positions and masses of point lenses */
double Re;			/* Einstein radius of the black in the formalism of Keeton */
double q, e;			/* axial ratio, x: major axis, y: minor axis */
double s;			/* core radius */
double q2, s2;			/* q^2, s^2 */
void findCritical();

/* image finding variables */
int nimages;
int nimagesMax;			/* maximum number of images */
double rmax=2.0;		/* maximum radius to search for images */
double phixx,  phixy,  phiyy;
double *xi, *yi;		/* image positions x & y */
double *mui;			/* image magnification */
double *abs_mui;			/* absoluteimage magnification */
void printLens();
void gslFindImage(), printImages();
void findImages();
double ran1(long *idum);
double xs, ys;			/* source positions */

long idum=-5;

int main()
{
  xl = (double *) malloc(nlenses * sizeof(double));
  yl = (double *) malloc(nlenses * sizeof(double));
  ml = (double *) malloc(nlenses * sizeof(double));


  xl[0] = -0.05; yl[0] = 0.0;
  xl[1] = +0.05; yl[1] = 0.0;
  Re = 0.0;
  ml[0] = Re*Re/2.0; ml[1] = ml[0]; /* equal mass */
  xs = 0.02;  ys = 0.04;

  s = 0.1;			/* core radius */
  q = 0.8;			/* axial ratio */
  q2 = q*q;
  s2 = s*s;
  e = sqrt(1.0-q*q);

  findCritical();
  findImages();
  printImages();

  return(0);
}

/* find all images */
double xs,  ys,  phix, phiy, mu;
void lensEquation(double x, double y);

/* obtain relative quantities in the lens equation at position (x, y)
 * source position (xs, ys),  phi_xx, phi_xy, phi_yy and the magnification.
 * All the above quantities are passed through a global array
 */
void lensEquation(double x, double y)
{
  int i;
  double r, r2;
  double dx, dy, ri2;
  double alphax, alphay;
  double Psi;
  double x2, y2;
  double denom;			/* denominator */

  x2 = x*x; y2 = y*y;
  r2 = x2+y2;  r = sqrt(r2);
  Psi = sqrt(q2*(s2+x2)+y2);

  alphax = 1.0/e *  atan(e*x/(Psi+s)) ;
  alphay = 1.0/e * atanh(e*y/(Psi+q2*s));

  phix = x-alphax;
  phiy = y-alphay;

  denom = 1.0/( Psi*( s2+q2*s2+r2+2.0*s*Psi ));
  phixx = 1.0-(q2*s2+y2+s*Psi) * denom;
  phixy = x*y * denom;
  phiyy = 1.0-(      s2+x2+s*Psi) * denom;

  mu = 1.0/(phixx*phiyy - phixy*phixy);
}


#define N 1000
void divide(double x, double y, double dx, double dy, int nlevels, FILE *fp);

void findCritical()
{
  int i, j;
  double mag[N][N];

  double xmin=-1.5, xmax=1.5, dx, x;
  double ymin=-1.5, ymax=1.5, dy, y;

  FILE *fp;

  dx = (xmax-xmin)/(double)(N-1);
  dy = (ymax-ymin)/(double)(N-1);

  for (i=0; i<N; i++) {
    x = xmin+dx*i;
    for (j=0; j<N; j++) {
      y = ymin+dy*j;

      lensEquation(x, y);
      mag[i][j] = mu;
    }
  }

  fp = fopen("cau.dat", "w");
  for (i=0; i<N-1; i++) {
    for (j=0; j<N-1; j++) {
      int nlevels=0;

      /* check parities of the four corners; if one or more are different, then subdivide the pixel  */
      if (mag[i][j]*mag[i][j+1]<0.0 || mag[i][j]*mag[i+1][j+1] < 0.0 || mag[i][j]*mag[i+1][j] <0.0) {
	divide(xmin+dx*i, ymin+dy*j, dx, dy, nlevels, fp);
      }
    }
  }
  fclose(fp);
}

/* find the critical curve recursively */
void divide(double x, double y, double dx, double dy, int nlevels, FILE *fp) {
  int i, j;
  double mag[3][3];

  if (nlevels >= 5) {
    lensEquation(x, y);
    fprintf(fp, "%f %f %f %f\n", x, y, phix, phiy);
    return;
  }

  /* find the magnifications for the nine points
   *    0 x 0
   *    x x x
   *    0 x 0
   * in principle, the four corner magnifications (marked as '0') are already known, and thus don't need to be re-computed. However, for coding simplicity, we recalculate their magnifications and those for the five points marked by 'x'
   */
  for (i=0; i<3; i++) {
    for (j=0; j<3; j++) {
      lensEquation(x+dx/2.0*i, y+dy/2.0*j);
      mag[i][j] = mu;
    }
  }

  /* check the parities of the four sub-squares */
  for (i=0; i<2; i++) {
    for (j=0; j<2; j++) {
      if (mag[i][j]*mag[i][j+1]<0.0 || mag[i][j]*mag[i+1][j+1] < 0.0 || mag[i][j]*mag[i+1][j] <0.0) {
	divide(x+dx/2.0*i, y+dy/2.0*j, dx*0.5, dy*0.5, nlevels+1, fp);
      }
    }
  }
}

#define EPS 1.0e-7
void findImages()
{
  size_t i;
  size_t j;
  double r, theta;
  double xi0, yi0;
  int nimages0=0;

  nimagesMax = 5*nlenses+1;	/* a guess of the maximum number of images */

  xi  = (double *) malloc(nimagesMax * sizeof(double));
  yi  = (double *) malloc(nimagesMax * sizeof(double));
  mui = (double *) malloc(nimagesMax * sizeof(double));
  abs_mui = (double *) malloc(nimagesMax * sizeof(double));

  /* search for the main positive parity image */
  xi0= xs/sqrt(xs*xs+ys*ys+TINY) * (1.0+sqrt(xs*xs+ys*ys));
  yi0= ys/sqrt(xs*xs+ys*ys+TINY) * (1.0+sqrt(xs*xs+ys*ys));

  gslFindImage(xi0, yi0);
  printf("GSL:initial + primary image: %f %f %d\n", xi0, yi0, nimages);

  /* search for the main negative parity image */
  xi0= xs/sqrt(xs*xs+ys*ys+TINY) * (-1.0+sqrt(xs*xs+ys*ys));
  yi0= ys/sqrt(xs*xs+ys*ys+TINY) * (-1.0+sqrt(xs*xs+ys*ys));

  gslFindImage(xi0, yi0);
  printf("GSL:initial + primary image: %f %f %d\n", xi0, yi0, nimages);

  for (i=0; i<nimages; i++) {
    printf("primary image: %f %f %f\n", xi[i], yi[i], mui[i]);
  }

  nimages0 = nimages;

  /* search for images close to the point lenses, we use 10 random
   * positions close to the lens to find the image
   */
  for (j=0; j<10; j++) {
    for (i=0; i<nlenses; i++) {
      xi0 = xl[i]+ml[i]*(2.0*ran1(&idum)-1.0);
      yi0 = yl[i]+ml[i]*(2.0*ran1(&idum)-1.0);

      gslFindImage(xi0, yi0);
    }
  }
  printf("point lens images found: %d\n", nimages-nimages0);
  nimages0 = nimages;

  /* search for images on a grid */
  for (i=0; i<N; i++) {
    r =  rmax * ran1(&idum);
    theta = 2.0*M_PI * ran1(&idum);

    gslFindImage(r*cos(theta), r*sin(theta));
  }
  printf("grid search images found: %d maxImage=%d\n",
	 nimages-nimages0, nimagesMax);
}

/* is a found image is new? If yes, returns 1; if not, returns 0 */
static int isNewImage(double ximage, double yimage)
{
  int i;
  double distance;

  for (i=0; i<nimages; i++) {
    distance = fabs(  (xi[i]-ximage)*(xi[i]-ximage) +
		      (yi[i]-yimage)*(yi[i]-yimage));

    if (distance < 2.0*EPS) {	/* image has already been found */
      return(0);
    }
  }

  return(1);
}

/* print image positions */
void printImages()
{
  int i;
  FILE *fp;

  fprintf(stderr, "number of images: %d\n", nimages);
  fp = fopen("image.dat", "w");
  for (i=0; i<nimages; i++) {
    fprintf(fp, "%d %f %f %e\n", i, xi[i], yi[i], mui[i]);
  }
  fclose(fp);
}

/* check image number */
void checkImageNumber()
{
  if ((nimages - nlenses)/2 !=  (nimages-nlenses)*0.5) { /* we assume non-isothermal sphere */
    printf("WARNING: image number is incorrect\n");
  }
}

/* print point lens positions */
void printLens()
{
  int i;
  FILE *fp;

  fprintf(stderr, "number of lenses: %d\n", nlenses);
  fp = fopen("lens.dat", "w");
  fprintf(fp, "%d\n", nlenses);
  for (i=0; i<nlenses; i++) {
    fprintf(fp, "%f %f\n", xl[i], yl[i]);
  }
  fclose(fp);
}

int print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));

  return(0);
}

struct rparams {double a; double b;};

int func(const gsl_vector * x,
	 void *params,
	 gsl_vector * f)
{
//  double a = ((struct rparams *) params)->a;
//  double b = ((struct rparams *) params)->b;

  const double xi = gsl_vector_get (x, 0);
  const double yi = gsl_vector_get (x, 1);

  lensEquation(xi, yi);

  gsl_vector_set (f, 0, xs-phix);
  gsl_vector_set (f, 1, ys-phiy);

  return GSL_SUCCESS;
}


/* find the image position using GNU science library routines */
void gslFindImage(double xi0, double yi0)
{
  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;

  int status;
  size_t iter = 0;

  const size_t n = 2;
  struct rparams p = {1.0, 100.0};

  double ximage, yimage;

  gsl_multiroot_function f = {&func, n, &p};

  gsl_vector *x = gsl_vector_alloc (n);

  /* initial guess image position */
  gsl_vector_set (x, 0, xi0);
  gsl_vector_set (x, 1, yi0);

  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, 2);

  gsl_multiroot_fsolver_set (s, &f, x);

  /*   print_state (iter, s); */

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);

    /*     print_state (iter, s);*/

    if (status)   /* check if solver is stuck */
      break;

    status =
      gsl_multiroot_test_residual (s->f, 1e-8);
  }
  while (status == GSL_CONTINUE && iter < 1000);

  if (status == GSL_SUCCESS) {
    ximage = gsl_vector_get (s->x, 0);
    yimage = gsl_vector_get (s->x, 1);

    if (isNewImage(ximage, yimage)) {

	xi[nimages] = ximage; yi[nimages] = yimage; mui[nimages] = mu;
	abs_mui[nimages] = fabs(mui[nimages]);

	nimages++;

	/* image number should be at least 2 && <= nimagesMax) */
	assert(nimages <= nimagesMax || nimages <= 3);
    }
  }

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
}
