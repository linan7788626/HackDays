int sign(float x);
float deg2rad(float pha);
void xy_rotate(float *x_in,float *y_in,int nx,int ny,float xcen,float ycen,float pha,float *x_out,float *y_out);
void gauss_2d(float *x1,float *x2, int nx1,int nx2,float *par,float *res);
void tophat_2d(float *x1,float *x2, int nx1,int nx2,float *par,float *res);
void lq_nie(float *x1,float *x2,int nx1,int nx2,float *lpar,float *alpha1,float *alpha2);
//void refine_critical();
void find_critical_curve(float *mu,int nx1,int nx2,float* res);
