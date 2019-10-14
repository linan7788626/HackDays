#define REAL double
#define h 0.71
#define Om0 0.27
#define Ol0 0.73
#define Otot Om0+Ol0
#define Ok0 0.0
//#define w 0.0-1.0
#define rho_crit0 2.78e11 //M_sun Mpc^-3 *h*h
#define rho_bar0 rho_crit0*Om0
#define sigma8 0.801 //wmap 7th
//----------------------------------------------------------------------------
#define apr 206269.43				//1/1^{''}
#define vc 2.9970e5				//km/s
#define M_G 4.3e-9				//(Mpc/h)^1 (Msun/h)^-1 (km/s)^2 
#define H0 100.0				//km/s/(Mpc/h)
#define pc 3.085677e16				//meter
#define kpc 3.085677e19				//meter
#define Mpc 3.085677e22				//meter
#define Msun 1.98892e30				//kg
#define yr 31536000.0/365.0*365.25		//second
#define Gyr yr*1e9				//second
//REAL M_G = 6.67259e-11;			//m^3/kg/s^2
//REAL 1Gly = 9.461e26cm = 9.461e24m = 9.461e21km; 
//REAL 1Mpc = 3.08568e24cm = 3.261566ly;
//----------------------------------------------------------------------------
REAL frand (void);
REAL sign(REAL x);
REAL ez(REAL x,void *params);
REAL ezf(REAL x);
REAL ez_integral(REAL x);
REAL tzf(REAL x);
REAL tz(REAL x,void *params);
REAL Hz(REAL x);
REAL scale_factor(REAL x);
REAL Dh();
REAL Th();
REAL Tl(REAL x);
REAL age_z(REAL x);
REAL Dc(REAL x);
REAL Dm(REAL x);
REAL Dl(REAL x);
REAL Dp(REAL x1,REAL x2);
REAL DistMod(REAL x);
REAL Omz(REAL x);
REAL Olz(REAL x);
REAL Okz(REAL x);
REAL rho_crit(REAL x);
REAL rho_bar(REAL x);
REAL Da(REAL x);
REAL Da2(REAL x1,REAL x2);
REAL dv(REAL z);
int retdz12(REAL x1,REAL x2,REAL *dz12);
REAL sigma_crit(REAL x1,REAL x2);
//----------------------------------------------------------------------------
struct f_params_2{REAL a,b;};
struct f_params_3{REAL a,b,c;};
//-------------------------------Def GSL fun----------------------------------
REAL sp_func(REAL x, void *params);
REAL dp_func(REAL x, void *params);
REAL tp_func(REAL x, void *params);
REAL qp_func(REAL x, void *params);
//-------------------------------1 parameter----------------------------------
int dif1p(REAL x,REAL (*fp)(REAL,void*),REAL * result,REAL * abserr);
int int_u1p(REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int int_sd1p(REAL a,REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int int_su1p(REAL a,REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int int_l1p(REAL a,REAL b,REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int f_root1p(REAL a,REAL b,REAL (*fp)(REAL, void*),REAL * result);
//-------------------------------2 parameters---------------------------------
int dif2p(REAL x,REAL alpha,REAL (*fp)(REAL,void*),REAL * result,REAL * abserr);
int int_u2p(REAL alpha,REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int int_sd2p(REAL alpha,REAL a,REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int int_su2p(REAL alpha,REAL a,REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int int_l2p(REAL alpha,REAL a,REAL b,REAL (*fp)(REAL, void*),REAL * result,REAL * error);
int f_root2p(REAL alpha,REAL a,REAL b,REAL (*fp)(REAL, void*),REAL * result);
//-------------------------------3 parameters---------------------------------
int dif3p(REAL x,REAL p1,REAL p2,REAL (*fp)(REAL,void*),
REAL * result,REAL * abserr);
int int_u3p(REAL p1,REAL p2,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int int_sd3p(REAL p1,REAL p2,REAL a,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int int_su3p(REAL p1,REAL p2,REAL a,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int int_l3p(REAL p1,REAL p2,REAL a,REAL b,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int f_root3p(REAL p1,REAL p2,REAL a,REAL b,REAL (*fp)(REAL, void*),REAL * result);
//-------------------------------4 parameters---------------------------------
int dif4p(REAL x,REAL p1,REAL p2,REAL p3,REAL (*fp)(REAL,void*),
REAL * result,REAL * abserr);
int int_u4p(REAL p1,REAL p2,REAL p3,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int int_sd4p(REAL p1,REAL p2,REAL p3,REAL a,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int int_su4p(REAL p1,REAL p2,REAL p3,REAL a,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int int_l4p(REAL p1,REAL p2,REAL p3,REAL a,REAL b,REAL (*fp)(REAL, void*),
REAL * result,REAL * error);
int f_root4p(REAL p1,REAL p2,REAL p3,REAL a,REAL b,REAL (*fp)(REAL, void*),
REAL * result);
