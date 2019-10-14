typedef struct{
    double redshift;
    int id;
} order;

//double *read_from_fits(char source_fits[256], long *nsx, long *nsy);
//void save_to_fits(double *lensed_map, int nlx, int nly, char *lensed_fits);
double array_max(double *array_in,int N);
double array_min(double *array_in,int N);
int array_max_int(int *array_in,int N);
int array_min_int(int *array_in,int N);
void masked_equal(int *array_in,int *array_out,int N,int value);
void saved_inside(double *array_in, double *array_out,int npixels,double value1,double value2);
int get_ngals(int *array_in,int npixels);
int aina(int *array, int N, int value);
void get_sub_z(char *filename, int *idlist0,int ngals,double *zlist,int *idlist);
int *get_sub_IDs(int *subfield, int npixels,int ngals);
void rotate2d_coor(double *x1,double *x2,double xc1,double xc2, double theta,double *x1_out,double *x2_out,int Np);
void rotate2d_vec(double x1,double x2,double xc1,double xc2, double theta,double *x1_out,double *x2_out);
double max2(double in1,double in2);
void crop_images(double * HUDF,double x0,double y0,int xlen,int ylen,double * sHUDF);
int compare(const void *a, const void *b);
void sort_zid(int *idlist,double *zlist,int ngals,int *idlist_sorted,double *zlist_sorted);
void get_redshift_bins(double zl,double *zlist_sorted, int ngals, int npb, double *zblist,int *nbins,int *nzles);
int *crop_sub_images(double *tm,double *ti,int xlen,int ylen,double x0,double y0,double *mHUDF,double *iHUDF,int *Ngals);
void cook_sources(int xlen,int ylen,int ngals,int npb,int nbins,int nzles,double *mHUDF, double *iHUDF, int *idlist_sorted, double *zHUDF);
double *cal_z_slices(double *tm,double *ti,int xlen,int ylen,double x0,double y0, double zl, int npb,int *Nbins,int *Nzles, int *Ngals);
void get_z_slices(double * zslices,int nbins, int xlen,int ylen,double *zblist,double *slices);
