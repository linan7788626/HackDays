#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <fitsio.h>
/*#include <omp.h>*/

#include "mycosmology.h"

/*//--------------------------------------------------------------------*/
/*void inverse_cic_omp(double *in_map, double *posy1, double *posy2, double ysc1, double ysc2,double dsi, int nsx, int nsy, int nlx, int nly, double *out_map) {*/

	/*int index;*/
	/*int i1,j1,i,j;*/
	/*double xb1,xb2;*/
	/*double ww1,ww2,ww3,ww4,wx,wy;*/

/*#pragma omp parallel num_threads(8)	\*/
	/*shared(in_map,nlx,nly,ysc1,ysc2,dsi,nsx,nsy,out_map) \*/
	/*private(index,i,j,i1,j1,xb1,xb2,wx,wy,ww1,ww2,ww3,ww4)*/
	/*{*/
	/*double *out_map_sp;*/
	/*out_map_sp = (double *)calloc(nlx*nly,sizeof(double));*/
	/*#pragma omp for schedule(dynamic,16)*/

	/*for(i=0;i<nlx;i++) for(j=0;j<nly;j++) {*/

		/*index = i*nlx+j;*/

		/*xb1 = (posy1[index]-ysc1)/dsi+(double)nsx/2.0-0.5;*/
		/*xb2 = (posy2[index]-ysc2)/dsi+(double)nsy/2.0-0.5;*/

		/*i1 = (int)xb1;*/
		/*j1 = (int)xb2;*/

		/*wx = 1.-(xb1-(double)(i1));*/
		/*wy = 1.-(xb2-(double)(j1));*/

		/*ww1 = wx*wy;*/
		/*ww2 = wx*(1.0-wy);*/
		/*ww3 = (1.0-wx)*wy;*/
		/*ww4 = (1.0-wx)*(1.0-wy);*/

		/*if (i1<0||i1>nsx-2||j1<0||j1>nsy-2) continue;*/

		/*out_map_sp[index] = ww1*in_map[i1*nsx+j1]*/
						  /*+ ww2*in_map[i1*nsx+j1+1]*/
						  /*+ ww3*in_map[(i1+1)*nsx+j1]*/
						  /*+ ww4*in_map[(i1+1)*nsx+j1+1];*/
	/*}*/

	/*#pragma omp critical*/
	/*{*/
		/*for(i=0;i<nlx;i++) for(j=0;j<nly;j++) {*/
			/*out_map[i*nly+j] += out_map_sp[i*nly+j];*/
		/*}*/
	/*}*/
	/*free(out_map_sp);*/
	/*}*/
/*}*/
//--------------------------------------------------------------------
void Interplation_on_source_plane(REAL *source_map, REAL *posy1, REAL *posy2, REAL ysc1, REAL ysc2,REAL dsi, int nsx, int nsy, int nlx, int nly, REAL *lensed_map) {

	int i1,j1,i,j;
	int index;
	REAL xb1,xb2;
	REAL ww1,ww2,ww3,ww4,wx,wy;

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++) {

		index = i*nly+j;

		xb1 = (posy1[index]-ysc1)/dsi+(REAL)nsx/2.0-0.5;
		xb2 = (posy2[index]-ysc2)/dsi+(REAL)nsy/2.0-0.5;

		i1 = (int)xb1;
		j1 = (int)xb2;

		wx = 1.-(xb1-(REAL)(i1));
		wy = 1.-(xb2-(REAL)(j1));

		ww1 = wx*wy;
		ww2 = wx*(1.0-wy);
		ww3 = (1.0-wx)*wy;
		ww4 = (1.0-wx)*(1.0-wy);

		if ((i1<0)||(i1>nsx-2)||(j1<0)||(j1>nsy-2)) {
			/*lensed_map[index] = 0.0;*/
			continue;
		}

		lensed_map[index] = ww1*source_map[i1*nsy+j1]
						  + ww2*source_map[i1*nsy+j1+1]
						  + ww3*source_map[(i1+1)*nsy+j1]
						  + ww4*source_map[(i1+1)*nsy+j1+1];
	}
}
//----------------------------------------------------------------------------------
void sfits_to_lfits(char *source_fits,REAL *posx1, REAL *posx2, REAL *alpha1, REAL *alpha2,REAL ysc1,REAL ysc2, REAL dsi, int nlx, int nly, char *lensed_fits) {

    fitsfile *fptr_in;
	int status,  nfound, anynull;
    long naxes[2],fpixel,npixels, i, j, index;

    float nullval;
    status = 0;

    fits_open_file(&fptr_in, source_fits, READONLY, &status);
    fits_read_keys_lng(fptr_in, "NAXIS", 1, 2, naxes, &nfound, &status);

	long nsx = naxes[0];
	long nsy = naxes[1];


    npixels  = nsx*nsy;
    fpixel   = 1;
    nullval  = 0;

	REAL *source_map = (REAL *)calloc(npixels,sizeof(REAL));
    fits_read_img(fptr_in, TDOUBLE, fpixel, npixels, &nullval,
                  source_map, &anynull, &status);

    fits_close_file(fptr_in, &status);

	REAL *posy1 = (REAL *)calloc(nlx*nly,sizeof(REAL));
	REAL *posy2 = (REAL *)calloc(nlx*nly,sizeof(REAL));

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++){
		index = i*nlx+j;
		posy1[index] = posx1[index]-alpha1[index];
		posy2[index] = posx2[index]-alpha2[index];
	}

	REAL *lensed_map = (REAL *)calloc(nlx*nly,sizeof(REAL));
	Interplation_on_source_plane(source_map,posy1,posy2,ysc1,ysc2,dsi,nsx,nsy,nlx,nly,lensed_map);
	//inverse_cic_omp(source_map,posy1,posy2,ysc1,ysc2,dsi,nsx,nsy,nlx,nly,lensed_map);
    fitsfile *fptr_out;
    int fpixel_out = fpixel;
    long bitpix_out =  DOUBLE_IMG;
    long naxis_out = 2;
    long npixels_out = nlx*nly;
    long naxes_out[naxis_out];
    naxes_out[0] = nlx;
    naxes_out[1] = nly;

    status = 0;

    remove(lensed_fits);
    fits_create_file(&fptr_out, lensed_fits, &status);
    fits_create_img(fptr_out,  bitpix_out, naxis_out, naxes_out, &status);
    fits_write_img(fptr_out, TDOUBLE, fpixel_out, npixels_out, lensed_map, &status);
    fits_close_file(fptr_out, &status);

	free(posy1);
	free(posy2);
}
//----------------------------------------------------------------------------------
void smap_to_lmap(REAL *source_map, int nsx, int nsy, REAL *posx1, REAL *posx2, REAL *alpha1, REAL *alpha2,REAL ysc1,REAL ysc2, REAL dsi, int nlx, int nly, REAL *lensed_map) {

	int i,j,index;

	REAL *posy1 = (REAL *)calloc(nlx*nly,sizeof(REAL));
	REAL *posy2 = (REAL *)calloc(nlx*nly,sizeof(REAL));

	for(i=0;i<nlx;i++) for(j=0;j<nly;j++){
		index = i*nlx+j;
		posy1[index] = posx1[index]-alpha1[index];
		posy2[index] = posx2[index]-alpha2[index];
	}

	Interplation_on_source_plane(source_map,posy1,posy2,ysc1,ysc2,dsi,nsx,nsy,nlx,nly,lensed_map);
	/*inverse_cic_omp(source_map,posy1,posy2,ysc1,ysc2,dsi,nsx,nsy,nlx,nly,lensed_map);*/

	free(posy1);
	free(posy2);
}
//----------------------------------------------------------------------------------
REAL *read_from_fits(char *source_fits, long *nsx, long *nsy) {

    fitsfile *fptr_in;
	int status,  nfound, anynull;
    long naxes[2],fpixel,npixels;

    float nullval;
    status = 0;

    fits_open_file(&fptr_in, source_fits, READONLY, &status);
    fits_read_keys_lng(fptr_in, "NAXIS", 1, 2, naxes, &nfound, &status);

	*nsx = naxes[0];
	*nsy = naxes[1];

    npixels  = naxes[0]*naxes[1];
    fpixel   = 1;
    nullval  = 0;

	REAL *source_map = (REAL *)calloc(npixels,sizeof(REAL));
    fits_read_img(fptr_in, TDOUBLE, fpixel, npixels, &nullval,
                  source_map, &anynull, &status);

    fits_close_file(fptr_in, &status);

	return source_map;
}
//----------------------------------------------------------------------------------
void save_to_fits(REAL *lensed_map, int nlx, int nly, char *lensed_fits) {

    fitsfile *fptr_out;
    int fpixel_out = 1;
    long bitpix_out =  DOUBLE_IMG;
    long naxis_out = 2;
    long npixels_out = nlx*nly;
    long naxes_out[naxis_out];
    naxes_out[0] = nlx;
    naxes_out[1] = nly;

    int status = 0;

    remove(lensed_fits);
    fits_create_file(&fptr_out, lensed_fits, &status);
    fits_create_img(fptr_out,  bitpix_out, naxis_out, naxes_out, &status);
    fits_write_img(fptr_out, TDOUBLE, fpixel_out, npixels_out, lensed_map, &status);
    fits_close_file(fptr_out, &status);
}
////----------------------------------------------------------------------------------
//void simple_sfits_to_lfits(char *source_fits,REAL *posx1, REAL *posx2, REAL *alpha1, REAL *alpha2,REAL ysc1,REAL ysc2, REAL dsi, int nlx, int nly, char *lensed_fits) {
//
//    fitsfile *fptr_in;
//	int status,  nfound, anynull;
//    long naxes[2],fpixel,npixels;
//
//    float nullval;
//    status = 0;
//
//    fits_open_file(&fptr_in, source_fits, READONLY, &status);
//    fits_read_keys_lng(fptr_in, "NAXIS", 1, 2, naxes, &nfound, &status);
//
//	long nsx = naxes[0];
//	long nsy = naxes[1];
//
//    npixels  = nsx*nsy;
//    fpixel   = 1;
//    nullval  = 0;
//
//	REAL *source_map = calloc(npixels,sizeof(REAL));
//    fits_read_img(fptr_in, TDOUBLE, fpixel, npixels, &nullval,
//                  source_map, &anynull, &status);
//
//    fits_close_file(fptr_in, &status);
//
//	REAL *lensed_map = calloc(nlx*nly,sizeof(REAL));
//	smap_to_lmap(source_map, nsx, nsy, posx1, posx2, alpha1, alpha2, ysc1, ysc2, dsi, nlx, nly, lensed_map);
//	save_to_fits(lensed_map, nlx, nly, lensed_fits);
//
//	free(source_map);
//	free(lensed_map);
//}
//----------------------------------------------------------------------------------
void simple_sfits_to_lfits(char *source_fits,REAL *posx1, REAL *posx2, REAL *alpha1, REAL *alpha2,REAL ysc1,REAL ysc2, REAL dsi, int nlx, int nly, char *lensed_fits) {

	long nsx = 0;
	long nsy = 0;
	REAL *source_map;
	source_map = read_from_fits(source_fits, &nsx, &nsy);

	REAL *lensed_map = (REAL *)calloc(nlx*nly,sizeof(REAL));
	smap_to_lmap(source_map, nsx, nsy, posx1, posx2, alpha1, alpha2, ysc1, ysc2, dsi, nlx, nly, lensed_map);
	save_to_fits(lensed_map, nlx, nly, lensed_fits);

	free(source_map);
	free(lensed_map);
}
//----------------------------------------------------------------------------------
