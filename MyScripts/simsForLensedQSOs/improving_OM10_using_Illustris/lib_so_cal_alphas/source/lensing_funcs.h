void kernel_green_iso(int, REAL *, REAL);
void kappa_to_phi(REAL *, REAL *, int, REAL);
void kernel_alphas_iso(int, REAL *, REAL *, REAL);
void kappa_to_alphas(REAL *, REAL *, REAL *, int, REAL);
void kernel_shears_iso(int, REAL *, REAL *, REAL);
void kappa_to_shears(REAL *, REAL *, REAL *, int, REAL);
void kernel_smooth(REAL,int, REAL *, REAL);
void smoothing_kappa(REAL *, REAL *, REAL *, int, REAL);
void kappa_to_kappac(REAL *, REAL *, int, REAL);
void nkappa_to_alphas(REAL *, int, REAL, REAL *, REAL *);
void calculate_mu(REAL *, REAL *, REAL *, int, REAL *);
void alphas_to_mu(REAL *,REAL *,int, REAL, REAL *);
//-------------------------------------------------------------
void sdens_to_alphas(REAL *, int, REAL, REAL, REAL, REAL *, REAL *);
void sdens_to_shears(REAL *, int, REAL, REAL, REAL, REAL *, REAL *);
void sdens_to_kappac(REAL *, int, REAL, REAL, REAL, REAL *);
void sdens_to_mu(REAL *, int, REAL, REAL, REAL, REAL *);
//-------------------------------------------------------------
void write_1_signal(char *, REAL *, int);
void print_1_signal(REAL *, int);
void write_2_signals(char *, char *, REAL *, REAL *, int);
void print_2_signals(REAL *, REAL *, int);
//-------------------------------------------------------------
void Loadin_sigma_mesh(int,char *,REAL *);
void Loadin_grids_mesh(REAL, REAL, REAL, int, REAL *, REAL *);
void sdens_to_kappa(REAL, REAL *, int, REAL, REAL, REAL, REAL *);
//-------------------------------------------------------------
void zero_padding(REAL *, int, int, REAL *out);
void print_matrix(REAL *, int, int);
void center_matrix(REAL *, int, int, REAL *);
