#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>

#include "mycosmology.h"
#include "ray_tracing_funcs.h"


typedef struct{
    double redshift;
    int id;
} order;

//double *read_from_fits(char source_fits[256], long *nsx, long *nsy) {
//
//	fitsfile *fptr_in;
//	int status,  nfound, anynull;
//	long naxes[2],fpixel,npixels;
//
//	float nullval;
//	status = 0;
//
//	fits_open_file(&fptr_in, source_fits, READONLY, &status);
//	fits_read_keys_lng(fptr_in, "NAXIS", 1, 2, naxes, &nfound, &status);
//
//	*nsx = naxes[0];
//	*nsy = naxes[1];
//
//	npixels  = naxes[0]*naxes[1];
//	fpixel	 = 1;
//	nullval  = 0;
//
//	double *source_map = calloc(npixels,sizeof(double));
//	fits_read_img(fptr_in, TDOUBLE, fpixel, npixels, &nullval,
//				  source_map, &anynull, &status);
//
//	fits_close_file(fptr_in, &status);
//
//	return source_map;
//}
////----------------------------------------------------------------------------------
//void save_to_fits(double *lensed_map, int nlx, int nly, char *lensed_fits) {
//
//	fitsfile *fptr_out;
//	int fpixel_out = 1;
//	long bitpix_out =  DOUBLE_IMG;
//	long naxis_out = 2;
//	long npixels_out = nlx*nly;
//	long naxes_out[naxis_out];
//	naxes_out[0] = (long)nlx;
//	naxes_out[1] = (long)nly;
//
//	int status = 0;
//
//	remove(lensed_fits);
//	fits_create_file(&fptr_out, lensed_fits, &status);
//	fits_create_img(fptr_out,  bitpix_out, naxis_out, naxes_out, &status);
//	fits_write_img(fptr_out, TDOUBLE, fpixel_out, npixels_out, lensed_map, &status);
//	fits_close_file(fptr_out, &status);
//}
//
////----------------------------------------------------------------------------------
double array_max(double *array_in,int N) {
	int i;
	double max=0;
	for (i=0;i<N;i++) {
		if (array_in[i]>max) {
			max = array_in[i];
		}
	}
	return max;
}

double array_min(double *array_in,int N) {
	int i;
	double min=0;
	for (i=0;i<N;i++) {
		if (array_in[i]<min) {
			min = array_in[i];
		}
	}
	return min;
}

int array_max_int(int *array_in,int N) {
	int i, max=-9999;
	for (i=0;i<N;i++) {
		if (max < array_in[i]) {
			max = array_in[i];
		}
	}
	return max;
}

int array_min_int(int *array_in,int N) {
	int i, min=0;
	for (i=0;i<N;i++) {
		if (array_in[i]<min) {
			min = array_in[i];
		}
	}
	return min;
}
//----------------------------------------------------------------------------------

void masked_equal(int *array_in,int *array_out,int N,int value) {
	int i;
	for (i=0;i<N;i++) {
		if (array_in[i] == value) {
			array_out[i] = 0;
		}
		else {
			array_out[i] = array_in[i];
		}
	}
}

void saved_inside(double *array_in, double *array_out,int npixels,double value1,double value2) {
	int i;
	for (i=0;i<npixels;i++) {
		if ((array_in[i]> value1) & (array_in[i] <value2)) {
			array_out[i] = 1.0;
		}
		else {
			array_out[i] = 0.0;
		}
	}
}

int get_ngals(int *array_in,int npixels) {
	int ngals = 0;

	int i;
	int *array_tmp = (int *)malloc(npixels*sizeof(int));
	for (i=0;i<npixels;i++) {
		array_tmp[i] = array_in[i];
	}

	int maxID = 0;
	maxID = array_max_int(array_tmp,npixels);
	while (maxID > 0) {
		masked_equal(array_tmp,array_tmp,npixels,maxID);
		ngals = ngals+1;
		maxID = array_max_int(array_tmp,npixels);
	}

	free(array_tmp);
	return ngals;
}

int aina(int *array, int N, int value) {
	int i;
	int res = 0;
	for (i=0;i<N;i++) {
		if (value == array[i]) {
			res = 1;
		}
	}
	return res;
}
//----------------------------------------------------------------------------------
void get_sub_z(char *filename, int *idlist0,int ngals,double *zlist,int *idlist) {

	FILE * file = fopen(filename, "r"); /* should check the result */
	char line[256];
	char *zstr_tmp;
	char *ind_tmp;
	//char *zstr_tmp = malloc(sizeof(char)*256);
	//char *ind_tmp;
	int id;

	int i=0;
	while (fgets(line, sizeof(line), file)) {
		ind_tmp = strtok(line,"\t");
		id = atoi(ind_tmp);
		zstr_tmp = strtok(NULL,"\t");
		if (aina(idlist0,ngals,id)) {
			while(ind_tmp !=NULL) {
				zstr_tmp=ind_tmp;
				ind_tmp = strtok(NULL,"\t");
			}
			zstr_tmp = strtok(zstr_tmp,"\n");
			zlist[i] = (double)atof(zstr_tmp);
			idlist[i] = id;
			i++;
		}
	}
	//free(zstr_tmp);
	fclose(file);
}


int *get_sub_IDs(int *subfield, int npixels,int ngals) {

	int *idlist = (int *)calloc(ngals,sizeof(int));

	int *tmp = (int *)malloc(npixels*sizeof(int));
	int i = 0;
	for (i=0;i<npixels;i++) {
		tmp[i] = subfield[i];
	}

	i = 0;
	int maxID = 0;
	maxID = array_max_int(tmp,npixels);
	for (i=0;i<ngals;i++) {
		masked_equal(tmp,tmp,npixels,maxID);
		idlist[i] = maxID;
		maxID = array_max_int(tmp,npixels);
	}

	free(tmp);
	return idlist;
}
////----------------------------------------------------------------------------------
//void rotate2d_coor(double *x1,double *x2,double xc1,double xc2, double theta,double *x1_out,double *x2_out,int Np) {
//	double Q[2][2];
//
//	Q[0][0] = cos(theta);
//	Q[0][1] = -1*sin(theta);
//	Q[1][0] = sin(theta);
//	Q[1][1] = cos(theta);
//
//	int i;
//	for (i=0; i < Np; i++) {
//		x1_out[i]  = (x1[i]-xc1)*Q[0][0] + (x2[i]-xc2)*Q[0][1];
//		x2_out[i]  = (x1[i]-xc1)*Q[1][0] + (x2[i]-xc2)*Q[1][1];
//	}
//}
////----------------------------------------------------------------------------------
void rotate2d_vec(double x1,double x2,double xc1,double xc2, double theta,double *x1_out,double *x2_out) {
	double Q[2][2];

	Q[0][0] = cos(theta);
	Q[0][1] = -1*sin(theta);
	Q[1][0] = sin(theta);
	Q[1][1] = cos(theta);

	*x1_out = (x1-xc1)*Q[0][0] + (x2-xc2)*Q[0][1];
	*x2_out = (x1-xc1)*Q[1][0] + (x2-xc2)*Q[1][1];
}

double max2(double in1,double in2) {
	if(in1 >= in2) {
		return in1;
	}
	else {
		return in2;
	}
}
//----------------------------------------------------------------------------------
void crop_images(double * HUDF,double x0,double y0,int xlen,int ylen,double * sHUDF) {

	double xco,yco;
	double imageAngle = 0.746333574;
	double imageWidth = 3363*2;
	double buff;
	buff = max2((double)xlen,(double)ylen)*sqrt(2.0) + 600;

	rotate2d_vec(x0-0.5,y0-0.5,0.0,0.0,imageAngle,&xco,&yco);

	int xl = 10500;
	int yl = 10500;

	int xci,yci;
	xci = (int)(xco*(imageWidth - buff)+xl/2);
	yci = (int)(yco*(imageWidth - buff)+yl/2);

	//printf("%d-------%d\n", xci,yci);

	//xci = 4787;
	//yci = 6652;

	int xhf = (int)(xlen/2.0);
	int yhf = (int)(ylen/2.0);

	int i,j;
	int index1,index2;
	for (i=0;i<(int)xlen;i++) for (j=0;j<(int)ylen;j++) {
		index1 = i*ylen+j;
		index2 = (yci+1-xhf+i)*yl+(xci+1-yhf+j);
		sHUDF[index1] = HUDF[index2];
	}
}

int compare(const void *a, const void *b) {
    //order *ia = (order *)a;
    //order *ib = (order *)b;
    //return (int)(1000.f*ia->redshift - 1000.f*ib->redshift);

	const order *ia = (order *)a;
	const order *ib = (order *)b;
	if (ia->redshift < ib->redshift)
	    return -1;
	else if (ia->redshift > ib->redshift)
	    return +1;
	else
	    return 0;

	/* float comparison: returns negative if b > a
	 *	and positive if a > b. We multiplied result by 100.0
	 *		to preserve decimal fraction */

}

void sort_zid(int *idlist,double *zlist,int ngals,int *idlist_sorted,double *zlist_sorted) {
	int i;

	order *list;
	list = (order *)malloc(ngals*sizeof(order));

	for (i=0; i<ngals; i++) {
		list[i].id = idlist[i];
		list[i].redshift = zlist[i];
	}

	qsort(list, ngals, sizeof(order), compare);

	for (i=0; i<ngals; i++) {
		zlist_sorted[i] = list[i].redshift;
		idlist_sorted[i] = list[i].id;
		//printf("-------------%lf----------------%d\n", zlist_sorted[i],idlist_sorted[i]);
	}
}

//void get_redshift_bins(double zl,double *zlist_sorted, int ngals, int npb, double *zblist,int *nbins,int *nzles) {
//
//	int Nzles = 0;
//	int i,j,index1=0,index2=0;
//	for (i=0;i<ngals;i++) {
//		if(zlist_sorted[i]<zl) {
//			Nzles++;
//		}
//	}
//
//	int Nbins;
//
//	if (npb%2==0) {
//		if (((ngals-Nzles)%npb) == 0) {
//			Nbins = (ngals-Nzles)/npb;
//			for (i=0;i<Nzles;i++) {
//				zblist[i] = 0;
//			}
//			for (i=0;i<Nbins;i++) {
//				for (j=0;j<npb;j++) {
//					index1 = Nzles+i*npb+j;
//					index2 = Nzles+i*npb+npb/2-1;
//					zblist[index1] = 0.5*(zlist_sorted[index2]+zlist_sorted[index2+1]);
//				}
//			}
//		}
//		else {
//			Nbins = (ngals-Nzles)/npb+1;
//			for (i=0;i<Nzles;i++) {
//				zblist[i] = 0;
//			}
//			for (i=0;i<Nbins-1;i++) {
//				for (j=0;j<npb;j++) {
//					index1 = Nzles+i*npb+j;
//					index2 = Nzles+i*npb+npb/2-1;
//					zblist[index1] = 0.5*(zlist_sorted[index2]+zlist_sorted[index2+1]);
//				}
//			}
//
//			for (j=(ngals-Nzles)%npb-1;j>=0;j--) {
//				index1 = ngals-1-j;
//				index2 = ngals-1-(ngals-npb*Nbins)/2;
//				zblist[index1] = 0.5*(zlist_sorted[index2]+zlist_sorted[index2+1]);
//			}
//		}
//	}
//	else {
//		if (((ngals-Nzles)%npb) == 0) {
//			Nbins = (ngals-Nzles)/npb;
//			for (i=0;i<Nzles;i++) {
//				zblist[i] = 0;
//			}
//			for (i=0;i<Nbins;i++) {
//				for (j=0;j<npb;j++) {
//					index1 = Nzles+i*npb+j;
//					index2 = Nzles+i*npb+(npb-1)/2;
//					zblist[index1] = zlist_sorted[index2];
//				}
//			}
//		}
//		else {
//			Nbins = (ngals-Nzles)/npb+1;
//			for (i=0;i<Nzles;i++) {
//				zblist[i] = 0;
//			}
//			for (i=0;i<Nbins-1;i++) {
//				for (j=0;j<npb;j++) {
//					index1 = Nzles+i*npb+j;
//					index2 = Nzles+i*npb+(npb-1)/2;
//					zblist[index1] = zlist_sorted[index2];
//				}
//			}
//
//			for (j=(ngals-Nzles)%npb-1;j>=0;j--) {
//				index1 = ngals-1-j;
//				if ((ngals-npb*Nbins)%2==0) {
//					index2 = Nzles+i*npb+npb/2-1;
//					index2 = ngals-1-(ngals-npb*Nbins)/2;
//					zblist[index1] = 0.5*(zlist_sorted[index2]+zlist_sorted[index2+1]);
//				}
//				else {
//					index2 = ngals-1-(ngals-npb*Nbins-1)/2;
//					zblist[index1] = zlist_sorted[index2];
//				}
//			}
//		}
//	}
//
//	*nbins=Nbins;
//	*nzles=Nzles;
//}

void get_redshift_bins(double zl,double *zlist_sorted, int ngals, int npb, double *zblist,int *nbins,int *nzles) {

	int Nzles = 0;
	int i,j,index1=0,index2=0;
	for (i=0;i<ngals;i++) {
		if(zlist_sorted[i]<zl) {
			Nzles++;
		}
	}

	int Nbins;

	if (((ngals-Nzles)%npb) == 0) {
		Nbins = (ngals-Nzles)/npb;
		for (i=0;i<Nzles;i++) {
			zblist[i] = 0;
		}
		for (i=0;i<Nbins;i++) {
			for (j=0;j<npb;j++) {
				index1 = Nzles+i*npb+j;
				if (npb%2==0) {
					index2 = Nzles+i*npb+npb/2-1;
					zblist[index1] = 0.5*(zlist_sorted[index2]+zlist_sorted[index2+1]);
				}
				else {
					index2 = Nzles+i*npb+(npb-1)/2;
					zblist[index1] = zlist_sorted[index2];
				}
			}
		}
	}
	else {
		int nrest;
		Nbins = (ngals-Nzles)/npb+1;
		for (i=0;i<Nzles;i++) {
			zblist[i] = 0;
		}
		for (i=0;i<Nbins-1;i++) {
			for (j=0;j<npb;j++) {
				index1 = Nzles+i*npb+j;
				if (npb%2==0) {
					index2 = Nzles+i*npb+npb/2-1;
					zblist[index1] = 0.5*(zlist_sorted[index2]+zlist_sorted[index2+1]);
				}
				else {
					index2 = Nzles+i*npb+(npb-1)/2;
					zblist[index1] = zlist_sorted[index2];
				}
			}
		}

		nrest = (ngals-Nzles)%npb;
		for (j=0;j<nrest;j++) {
			index1 = Nzles+(Nbins-1)*npb+j;
			if ((nrest)%2==0) {
				index2 = Nzles+(Nbins-1)*npb+nrest/2-1;
				zblist[index1] = 0.5*(zlist_sorted[index2]+zlist_sorted[index2+1]);
			}
			else {
				index2 = Nzles+(Nbins-1)*npb+(nrest-1)/2;
				zblist[index1] = zlist_sorted[index2];
			}
		}
	}
	*nbins=Nbins;
	*nzles=Nzles;
}

int *crop_sub_images(double *tm,double *ti,int xlen,int ylen,double x0,double y0,double *mHUDF,double *iHUDF,int *Ngals) {

	crop_images(tm,x0,y0,xlen,ylen,mHUDF);
	crop_images(ti,x0,y0,xlen,ylen,iHUDF);

	int *mm;
	int i;
	mm = (int *)calloc(xlen*ylen,sizeof(int));
	for (i=0;i<xlen*ylen;i++) {
		mm[i] = (int)(mHUDF[i]);
	}

	int ngals = get_ngals(mm,xlen*ylen);

	*Ngals = ngals;

	int *idlist0;
	idlist0 = get_sub_IDs(mm,xlen*ylen,ngals);

	free(mm);

	return idlist0;
}

void cook_sources(int xlen,int ylen,int ngals,int npb,int nbins,int nzles,double *mHUDF, double *iHUDF, int *idlist_sorted, double *zHUDF) {


	double * masks = (double *)calloc(xlen*ylen,sizeof(double));
	double * masks_tmp = (double *)calloc(xlen*ylen,sizeof(double));

	int i,j,k,l,index1=0,index2=0;
	for (l=0;l<nzles;l++) {
		saved_inside(mHUDF,masks_tmp,xlen*ylen,idlist_sorted[l]-0.1,idlist_sorted[l]+0.1);
		for(i=0;i<xlen;i++) for(j=0;j<ylen;j++) {
			index1 = i*ylen+j;
			masks[index1] = masks[index1] + masks_tmp[index1];
		}
	}

	for(i=0;i<xlen;i++) for(j=0;j<ylen;j++) {
		index1 = i*ylen+j;
		zHUDF[index1] = zHUDF[index1] + iHUDF[index1]*masks[index1];
	}

	//char fits_out[256];
	//double *source_tmp = (double *)calloc(xlen*ylen,sizeof(double));
	//for(j=0;j<xlen*ylen;j++) {
	//	source_tmp[j]=source_tmp[j]+zHUDF[j];
	//	//source_tmp[j]=source_tmp[j]+masks_tmp[j];
	//	//source_tmp[j]=source_tmp[j]+iHUDF[j];
	//}
	//sprintf(fits_out, "./test/sources.fits");
	//save_to_fits(source_tmp,xlen,ylen,fits_out);
	//free(source_tmp);

	free(masks);
	free(masks_tmp);

	for (k=0;k<nbins;k++) {
		double * masks = (double *)calloc(xlen*ylen,sizeof(double));
		double * masks_tmp = (double *)calloc(xlen*ylen,sizeof(double));

		for (l=0;l<npb;l++) {
			if ((nzles+k*npb+l)>ngals-1) continue;
			saved_inside(mHUDF,masks_tmp,xlen*ylen,idlist_sorted[nzles+k*npb+l]-0.1,idlist_sorted[nzles+k*npb+l]+0.1);
			for(i=0;i<xlen;i++) for(j=0;j<ylen;j++) {
				index2 = i*ylen+j;
				masks[index2] = masks[index2] + masks_tmp[index2];
			}
		}

		for(i=0;i<xlen;i++) for(j=0;j<ylen;j++) {
			index1 = (k+1)*xlen*ylen+(i*ylen+j);
			index2 = i*ylen+j;
			zHUDF[index1] = zHUDF[index1] + iHUDF[index2]*masks[index2];
		}

		//char fits_out[256];
		//double *source_tmp = (double *)calloc(xlen*ylen,sizeof(double));
		//for(j=0;j<xlen*ylen;j++) {
		//	source_tmp[j]=source_tmp[j]+zHUDF[(1)*xlen*ylen+j];
		//	source_tmp[j]=source_tmp[j]+masks_tmp[j];
		//	//source_tmp[j]=source_tmp[j]+iHUDF[j];
		//}
		//sprintf(fits_out, "./test/sources.fits");
		//save_to_fits(source_tmp,xlen,ylen,fits_out);
		//free(source_tmp);


		free(masks);
		free(masks_tmp);
	}

}
//----------------------------------------------------------------------------------
double *cal_z_slices(double *tm,double *ti,int xlen,int ylen,double x0,double y0, double zl, int npb,int *Nbins,int *Nzles, int *Ngals) {

	char * src_rsf = "/media/star2/nanli/piecewise_lensing_pipeline/sources/completeRedshifts.cat";

	int ngals = 0;
	int nbins = 0;
	int nzles = 0;

	double *mHUDF = (double *)malloc(xlen*ylen*sizeof(double));
	double *iHUDF = (double *)malloc(xlen*ylen*sizeof(double));

	int *idlist0;
	idlist0 =  crop_sub_images(tm,ti,xlen,ylen,x0,y0,mHUDF,iHUDF,&ngals);
	int *idlist = (int *)calloc(ngals,sizeof(int));
	double *zlist = (double *)calloc(ngals,sizeof(double));
	get_sub_z(src_rsf, idlist0, ngals, zlist,idlist);
	int *idlist_sorted = (int *)malloc(ngals*sizeof(int));
	double *zlist_sorted = (double *)malloc(ngals*sizeof(double));
	sort_zid(idlist,zlist,ngals,idlist_sorted,zlist_sorted);
	double *zblist = (double *)malloc(ngals*sizeof(double));
	get_redshift_bins(zl,zlist_sorted,ngals,npb,zblist,&nbins,&nzles);

	double *zHUDF = (double *)malloc((nbins+1)*xlen*ylen*sizeof(double));
	cook_sources(xlen,ylen,ngals,npb,nbins,nzles,mHUDF,iHUDF,idlist_sorted,zHUDF);

	double *z_slices = (double *)malloc(((nbins+1)*xlen*ylen+(nbins+1))*sizeof(double));
	int i;
	z_slices[0] = 0.0;
	for (i=1;i<(nbins+1);i++) {
		z_slices[i] = zblist[nzles+(i-1)*npb];
	}
	for (i=(nbins+1);i<((nbins+1)*xlen*ylen+(nbins+1));i++) {
		z_slices[i] = zHUDF[i-(nbins+1)];
	}

	//char fits_out[256];
	//double *source_tmp = (double *)calloc(xlen*ylen,sizeof(double));
	//for(i=0;i<xlen*ylen;i++) {
	//	source_tmp[i]=source_tmp[i]+zHUDF[i];
	//	//source_tmp[i]=source_tmp[i]+z_slices[i];
	//	//source_tmp[i]=source_tmp[i]+masks_tmp[i];
	//	//source_tmp[i]=source_tmp[i]+iHUDF[i];
	//}
	//sprintf(fits_out, "./test/sources.fits");
	//save_to_fits(source_tmp,xlen,ylen,fits_out);
	//free(source_tmp);

	*Ngals = ngals;
	*Nbins = nbins;
	*Nzles = nzles;

	free(idlist0);
	free(idlist);
	free(zlist);
	free(idlist_sorted);
	free(zlist_sorted);
	free(zblist);
	free(mHUDF);
	free(iHUDF);
	free(zHUDF);
	return z_slices;
}
//----------------------------------------------------------------------------------
void get_z_slices(double * zslices,int nbins, int xlen,int ylen,double *zblist,double *slices) {
	int i,j;
	for(i=0;i<nbins+1;i++) {
		zblist[i]= zslices[i];
		for(j=0;j<xlen*ylen;j++) {
			slices[i*xlen*ylen+j]= zslices[(nbins+1)+i*xlen*ylen+j];
		}
	}

	//for(i=0;i<nbins+1;i++) {
	//	for(j=0;j<xlen*ylen;j++) {
	//		slices[i*xlen*ylen+j]= zslices[(nbins+1)+i*xlen*ylen+j];
	//	}
	//}
}
//----------------------------------------------------------------------------------
//int main(int argc, char **argv) {
//
//
//    char * sourceM = argv[1];
//    char * sourceI = argv[1];
//
//	//int *idlist0;
//	int xlen = 1000;
//	int ylen = 1000;
//	double x0 = 0.562;
//	double y0 = 0.562;
//	double zl = 0.165;
//	int npb = 10;
//
//	int nbins;
//	int nzles;
//	int ngals;
//	double * z_slices;
//	z_slices = cal_z_slices(sourceM,sourceI,xlen,ylen,x0,y0,zl,npb,&nbins,&nzles,&ngals);
//
//	//printf("%d %d %d  -------------------------------\n",nbins,nzles,ngals);
//
//	//int i,j;
//    //char fits_out[256];
//	//double *source_tmp = (double *)calloc(xlen*ylen,sizeof(double));
//	//for(i=0;i<nbins+1;i++) {
//	//	for(j=0;j<xlen*ylen;j++) {
//	//		source_tmp[j]=source_tmp[j]+z_slices[(nbins+1)+i*xlen*ylen+j];
//	//	}
//	//	sprintf(fits_out, "test/sources.fits");
//	//	save_to_fits(source_tmp,xlen,ylen,fits_out);
//	//}
//
//    //char fits_out_tmp[256];
//	//for(i=0;i<nbins+1;i++) {
//	//	for(j=0;j<xlen*ylen;j++) {
//	//		source_tmp[j]= z_slices[(nbins+1)+i*xlen*ylen+j];
//	//	}
//	//	printf("%d %lf----------------//\n",i,z_slices[i]);
//	//	sprintf(fits_out_tmp, "test/source_%d.fits",i);
//	//	save_to_fits(source_tmp,xlen,ylen,fits_out_tmp);
//	//}
//
//	free(z_slices);
//	return 0;
//}

