typedef struct{
	float x;
	float y;
	float z;
} PARTICLE;

void Loadin_particle_main(long Np, char *fname, PARTICLE *particle);
void Loadin_particle_main_ascii(long Np, char *fname, PARTICLE *particle);
int findHsml(PARTICLE *particles, long *NumP, long * NgbP, double *BS, float * SmoothLength);

float si_weight(float x);
float  *sp_make_sph(long Nc,float bsz,float x,float y,float l,long *oi_l,long *oj_l,long *onbx,long *onby);
void pin_matrix(long Nc,long i_l,long j_l,long nbx,long nby,float *in1,float *in2,float *out);
void Make_cell_SPH(long Nc,float bsz,long Np, PARTICLE *particle, float * SmoothLength, float *sdens);

int cal_sph_sdens(char *in_part,float bsz,long  Nc,float dsx,long Ngb,long Np,float xc1,float xc2,float xc3,float mass_particle, float * posx1, float * posx2, float *sdens);
void write_3_signals(char *out1o,char *out2o,char *out3o,float *in1,float *in2, float *in3,int Nc, int ind);
