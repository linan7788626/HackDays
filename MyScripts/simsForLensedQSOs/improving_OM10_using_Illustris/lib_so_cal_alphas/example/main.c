#include <assert.h>
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <stdint.h>
#include <tiffio.h>

#define REAL double
/************* LENSING HEADERS ***************/
#include "mycosmology.h"
#include "lensing_funcs.h"
#include "ray_tracing_funcs.h"

//----------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    //*******************************************************************************
    // input parameters
    char *input_file = argv[1];
    int Nc = atoi(argv[2]);
    REAL boxsize = (REAL)atof(argv[3]);
    REAL zl = (REAL)atof(argv[4]);
    REAL zs = (REAL)atof(argv[5]);
	//REAL pmass = 1.85e9;
	//REAL pmass = 1.85e8;
	REAL pmass = 1.0;

	//// output files
    //char *out_posx1  = argv[6];
    //char *out_posx2  = argv[7];
    //char *out_alpha1 = argv[8];
    //char *out_alpha2 = argv[9];
    //char *out_shear1 = argv[10];
    //char *out_shear2 = argv[11];
    //char *out_kappa = argv[12];
    //char *out_mu = argv[13];

	//// output files
	//char *source_fits = argv[14];
	//char *lensed_fits = argv[15];
    ////*******************************************************************************
	//// Calculate all lensing signals
    ////*******************************************************************************
	//// Calculate sdens
	//REAL *sdens = calloc(Nc*Nc,sizeof(REAL));
	//Loadin_sigma_mesh(Nc*Nc, input_file, sdens);
    ////*******************************************************************************
	//// Call sdens_to_alphas
	//REAL *alpha1 = calloc(Nc*Nc,sizeof(REAL));
	//REAL *alpha2 = calloc(Nc*Nc,sizeof(REAL));
	//sdens_to_alphas(pmass, sdens, Nc, boxsize, zl, zs, alpha1, alpha2);
	//write_2_signals(out_alpha1,out_alpha2,alpha1,alpha2,Nc);
	////print_2_signals(alpha1,alpha2,Nc);
	//free(alpha1);
	//free(alpha2);
    ////*******************************************************************************
	//// Call sdens_to_shear
	//REAL *shear1 = calloc(Nc*Nc,sizeof(REAL));
	//REAL *shear2 = calloc(Nc*Nc,sizeof(REAL));
	//sdens_to_shears(pmass, sdens, Nc, boxsize, zl, zs, shear1, shear2);
	//write_2_signals(out_shear1,out_shear2,shear1,shear2,Nc);
	//free(shear1);
	//free(shear2);
    ////*******************************************************************************
	//// Call sdens_to_kappac
	//REAL *kappac = calloc(Nc*Nc,sizeof(REAL));
	//sdens_to_kappac(pmass, sdens, Nc, boxsize, zl, zs, kappac);
	//write_1_signal(out_kappa,kappac,Nc);
	//free(kappac);
    ////*******************************************************************************
	//// Calculate mu
	//REAL *mu = calloc(Nc*Nc,sizeof(REAL));
	//sdens_to_mu(pmass, sdens, Nc, boxsize, zl, zs, mu);
	//write_1_signal(out_mu,mu,Nc);
    ////*******************************************************************************
	//// Loadin posxs
	//REAL *posx1 = calloc(Nc*Nc,sizeof(REAL));
	//REAL *posx2 = calloc(Nc*Nc,sizeof(REAL));
	//Loadin_grids_mesh(boxsize, 0.0, 0.0, Nc, posx1, posx2);
	//write_2_signals(out_posx1,out_posx2,posx1,posx2,Nc);
	//free(posx1);
	//free(posx2);
    ////*******************************************************************************
	//// Free memory
	//free(sdens);



    ////*******************************************************************************
	//// Ray Tracing
    ////*******************************************************************************
	//// Calculate sdens
	//REAL *sdens = calloc(Nc*Nc,sizeof(REAL));
	//Loadin_sigma_mesh(Nc*Nc, input_file, sdens);
    ////*******************************************************************************
	//// Call sdens_to_alphas
	//REAL *alpha1 = calloc(Nc*Nc,sizeof(REAL));
	//REAL *alpha2 = calloc(Nc*Nc,sizeof(REAL));
	//sdens_to_alphas(pmass, sdens, Nc, boxsize, zl, zs, alpha1, alpha2);
    ////*******************************************************************************
	//// Loadin posxs
	//REAL *posx1 = calloc(Nc*Nc,sizeof(REAL));
	//REAL *posx2 = calloc(Nc*Nc,sizeof(REAL));
	//Loadin_grids_mesh(boxsize, 0.0, 0.0, Nc, posx1, posx2);

	//REAL ysc1 = -4.20;
	//REAL ysc2 = 3.10;
	//REAL dsi = 0.003;

	////simple_sfits_to_lfits(source_fits, posx1, posx2, alpha1, alpha2, ysc1, ysc2, dsi, Nc, Nc, lensed_fits);
	//sfits_to_lfits(source_fits, posx1, posx2, alpha1, alpha2, ysc1, ysc2, dsi, Nc, Nc, lensed_fits);

	//free(posx1);
	//free(posx2);
	//free(alpha1);
	//free(alpha2);
	//free(sdens);

    //*******************************************************************************
	// Combine with Python
    //*******************************************************************************
	// Calculate sdens
	REAL *sdens = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	Loadin_sigma_mesh(Nc*Nc, input_file, sdens);
    //*******************************************************************************
	// Call sdens_to_alphas
	REAL *alpha1 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	REAL *alpha2 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	sdens_to_alphas(pmass, sdens, Nc, boxsize, zl, zs, alpha1, alpha2);
	print_2_signals(alpha1,alpha2,Nc);

	free(alpha1);
	free(alpha2);
	free(sdens);
    return 0;
}
