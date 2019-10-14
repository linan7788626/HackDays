#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mycosmology.h"
#include "ray_tracing_funcs.h"
#include "cook_sources.h"
//----------------------------------------------------------------------------------
void Loadin_grids_mesh(REAL boxsize, REAL xc1, REAL xc2, int Ncc, REAL *posx1, REAL *posx2) {
    REAL dsx = boxsize/(REAL)Ncc;

    int i,j;
    int index;
    for (i=0; i<Ncc; i++) for (j=0; j<Ncc; j++) {
        index = i*Ncc+j;
        posx1[index] = dsx*(REAL)(i)-boxsize*0.5+0.5*dsx+xc1;
        posx2[index] = dsx*(REAL)(j)-boxsize*0.5+0.5*dsx+xc2;
    }
}

//----------------------------------------------------------------------------------
int ray_tracing_single(REAL *alpha1,REAL *alpha2,int Nc,REAL bsz,REAL zl,REAL *lensed_imgs_i) {

	double zs0=10.0;
	int xlen = 1024;
	int ylen = 1024;
	int npb = 10;
	REAL xsc1 = 0.8177;
	REAL xsc2 = 0.6234;
	double dsi = 0.03;

	char *sources_dir = "/media/star2/nanli/piecewise_lensing_pipeline/sources";
	char file_mask[256];
	sprintf(file_mask,"%s/h_udf_masks.fits",sources_dir);
  	char file_image_i[256];
	sprintf(file_image_i,"%s/h_udf_wfc_i_drz_img.fits",sources_dir);

	long nsx;
	long nsy;

	double * sourceM;
	double * sourceI;

	sourceM = read_from_fits(file_mask, &nsx, &nsy);
	sourceI = read_from_fits(file_image_i, &nsx, &nsy);

	double *posx1 = (double *)malloc(Nc*Nc*sizeof(double));
	double *posx2 = (double *)malloc(Nc*Nc*sizeof(double));
	Loadin_grids_mesh(bsz,0.0,0.0,Nc,posx1,posx2);
	//-------------------------------------------------------------------------------------
	int nbins=0;
	int nzles=0;
	int ngals=0;

	double *zslices_i;
	zslices_i = cal_z_slices(sourceM,sourceI,xlen,ylen,xsc1,xsc2,zl,npb,&nbins,&nzles,&ngals);

	double *zblist = (double *)malloc((nbins+1)*sizeof(double));

	double *slices_i = (double *)malloc((nbins+1)*xlen*ylen*sizeof(double));
	get_z_slices(zslices_i,nbins,xlen,ylen,zblist,slices_i);

	int i,k;
	double zs = 0;
	double factor_z = 0,factor_s = 0;
	double factor_z0 = Da(zs0)/Da2(zl,zs0);
	for (i=0;i<Nc*Nc;i++) {
		lensed_imgs_i[i] = 0.0;
	}

	for (k=0;k<nbins+1;k++) {
		double *source_tmp_i = (double *)malloc(xlen*ylen*sizeof(double));
		double *alpha1_tmp = (double *)malloc(Nc*Nc*sizeof(double));
		double *alpha2_tmp = (double *)malloc(Nc*Nc*sizeof(double));
		double *lensed_tmp_i = (double *)calloc(Nc*Nc,sizeof(double));

		zs = zblist[k];

		factor_z = zs>zl? Da2(zl,zs)/Da(zs)*factor_z0:0;
		factor_s = 1.0;

		for (i=0;i<xlen*ylen;i++) {
			source_tmp_i[i] = slices_i[xlen*ylen*(k)+i]*factor_s;
		}

		/*char fits_out_source_i[256];*/
		/*sprintf(fits_out_source_i, "/media/star2/nanli/piecewise_lensing_pipeline/slices/illustris_source_planes_%d_%lf_I.fits",k,zs);*/
		/*save_to_fits(source_tmp_i,xlen,ylen,fits_out_source_i);*/

		for (i=0;i<Nc*Nc;i++) {
			alpha1_tmp[i] = alpha1[i]*factor_z;
			alpha2_tmp[i] = alpha2[i]*factor_z;
		}

		smap_to_lmap(source_tmp_i,xlen,ylen,posx1,posx2,alpha1_tmp,alpha2_tmp,0.0,0.0,dsi,Nc,Nc,lensed_tmp_i);

		for (i=0;i<Nc*Nc;i++) {
			lensed_imgs_i[i] = lensed_imgs_i[i]+lensed_tmp_i[i];
		}

		/*sprintf(fits_out_source_i, "/media/star2/nanli/piecewise_lensing_pipeline/slices/lensed_illustris_source_planes_%d_%lf_I.fits",k,zs);*/
		/*save_to_fits(lensed_tmp_i,Nc,Nc,fits_out_source_i);*/

		free(source_tmp_i);
		free(alpha1_tmp);
		free(alpha2_tmp);
		free(lensed_tmp_i);
	}
	free(zslices_i);
	free(zblist);
	free(slices_i);
	free(sourceM);
	free(sourceI);
	free(posx1);
	free(posx2);
    return 0;
}
