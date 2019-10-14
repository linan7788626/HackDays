#define REAL double
//----------------------------------------------------------------------------------
void Loadin_grids_mesh(REAL boxsize, REAL xc1, REAL xc2, int Ncc, REAL *posx1, REAL *posx2);
//----------------------------------------------------------------------------------
void ms_ray_tracing(int Nc,double *posx1,double *posx2,double *alpha1,double *alpha2,int xlen,int ylen,double x0,double y0,double zl,int npb,double zs0,double dsi,double *sourceM,double *sourceI,double *sourceV,double *sourceB,double *lensed_images_i,double *lensed_images_v,double *lensed_images_b);
void stack_images(double *img1,int nx1,int ny1,double *img2,int nx2,int ny2,double *img_out);
//----------------------------------------------------------------------------------
void particles_to_file(char *out,PARTICLE *particle, int Np);
void write_3_signals(char *out1o,char *out2o,char *out3o,float *in1,float *in2, float *in3,int Nc, int ind);
int cal_sph_sdens(char *in_part,float bsz,long  Nc,REAL zl,long Ngb,long Np,float xc1,float xc2,float xc3,float mass_particle, REAL *sdens);
int produce_all_lensing(REAL * sdens, int Nc, REAL boxsize, REAL zl, REAL zs0, char * file_mask, char * file_image_i, char * file_image_v, char * file_image_b, char * gals_list, int i, REAL pmass);
int mpi_main(char *input1,char *input2,int Nc,REAL boxsize,REAL zl,REAL zs0,char * file_mask,char * file_image_i,char * file_image_v,char * file_image_b,char * lensed_dir,int ind);
int mpi_main2(REAL *alpha1,REAL *alpha2,int Nc,REAL boxsize,REAL zl,REAL zs0,char * file_mask,char * file_image_i,char * file_image_v,char * file_image_b,char * lensed_dir,int ind);
int mpi_main3(REAL *alpha1,REAL *alpha2,int Nc,REAL boxsize,REAL zl,REAL zs0,int xlen,int ylen,int npb,REAL xsc1,REAL xsc2,char * file_mask,char * file_image_i,char * file_image_v,char * file_image_b,char * lensed_dir,int ind);
int mpi_main_single(REAL *alpha1,REAL *alpha2,int Nc,REAL bsz,REAL zl,REAL *lensed_imgs_i);
