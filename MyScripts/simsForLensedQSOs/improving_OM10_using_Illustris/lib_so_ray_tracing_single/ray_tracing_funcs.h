void inverse_cic_omp(double *in_map, double *posy1, double *posy2, double ysc1, double ysc2,double dsi, int nsx, int nsy, int nlx, int nly, double *out_map);
void Interplation_on_source_plane(REAL *, REAL *, REAL *, REAL, REAL,REAL, int, int, int, int, REAL *);
void sfits_to_lfits(char *,REAL *, REAL *, REAL *, REAL *,REAL,REAL, REAL, int, int, char *);
void smap_to_lmap(REAL *, int, int, REAL *, REAL *, REAL *, REAL *, REAL, REAL, REAL, int, int, REAL *);
REAL *read_from_fits(char *, long *, long *);
void save_to_fits(REAL *, int, int, char *);
void simple_sfits_to_lfits(char *, REAL *, REAL *, REAL *, REAL *, REAL, REAL, REAL, int, int, char *);
