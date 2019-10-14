#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
//#include "non_lensing_funcs.h"
#include "mycosmology.h"
#include "fft_convolve.h"
//#include "constants.h"
//#include "fitsio.h"
//----------------------------------------------------------------------------------
void Loadin_sigma_mesh(int Num,char *fname,REAL* sdens) {
	int i;
	REAL x;
	//REAL sdens_max = 0.0;
	FILE *f_rb;

	//printf("load in density mesh from file...");
	f_rb = fopen(fname,"rb");

	int unused __attribute__((unused));
	for(i=0;i<Num;i++) {
		unused=fread(&x,sizeof(REAL),1,f_rb);
		sdens[i] = x;
		//sdens_max = (sdens[i]>sdens_max)? sdens[i]:sdens_max;
		//sdens_max = (sdens[i]<sdens_max)? sdens[i]:sdens_max;
	}
	//printf("Maxium sdens = %5.2f\n", sdens_max);
	//printf("%d %s\n", Num,fname);
	fclose(f_rb);
	//printf("done!\n");

}

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
void sdens_to_kappa(REAL* sdens, int Ncc, REAL Dcell, REAL zl, REAL zs, REAL *kappa) {

	REAL scrit;
	scrit = sigma_crit(zl,zs);//*(apr*apr/(Da(zl)*Da(zl)));

	int i,j,index;

	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		kappa[index] = sdens[index]/scrit;
	}
}
//--------------------------------------------------------------------
void zero_padding(REAL *in, int nx, int ny, REAL *out) {

	int i,j,index1,index2;
	for ( i = 0; i < nx; i++ ) {
		for ( j = 0; j < ny; j++ ) {
			index1 = i*ny+j;
			index2 = i*2*ny+j;
			out[index2] = in[index1];
		}
	}
}
//--------------------------------------------------------------------
void print_matrix(REAL *in, int nx, int ny) {

	int i,j,index;
	for ( i = 0; i < nx; i++ )
	{
		for ( j = 0; j < ny; j++ )
		{
			index = i*ny+j;
			printf ( "[%d,%d]%06f ",i, j, in[index]);
		}
		printf("\n");
	}
	printf("------------------------------------------\n");
}
//--------------------------------------------------------------------
void center_matrix(REAL *in, int nx, int ny, REAL *out) {

	int i,j,index1,index2;
	for (i=nx/4;i<3*nx/4; i++ ) {
		for (j=ny/4;j<3*ny/4; j++ ) {
			index1 = i*ny+j;
			index2 = (i-nx/4)*ny/2+(j-ny/4);
			out[index2] = in[index1];
		}
	}
}
//--------------------------------------------------------------------
void corner_matrix(REAL *in, int nx, int ny, REAL *out) {

	int i,j,index1,index2;
	for (i=0;i<nx/2; i++ ) {
		for (j=0;j<ny/2; j++ ) {
			index1 = i*ny+j;
			index2 = i*ny/2+j;
			out[index2] = in[index1];
		}
	}
}
//--------------------------------------------------------------------
void kernel_green_iso(int Ncc, REAL *in, REAL Dcell) {
	int i,j;
	REAL x,y,r;
	REAL epsilon = 0.00000001*Dcell;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <=(Ncc/2)  && j <=(Ncc/2)) {
			x = (REAL)(i)*Dcell-0.5*Dcell;
			y = (REAL)(j)*Dcell-0.5*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);

			if(r > Dcell*(REAL)Ncc/2.0) {
				in[i*Ncc+j] = 0.0;
			}
			else {
				in[i*Ncc+j] = 1.0/M_PI*log(r);
			}
		}
		else {
			if(i <= Ncc/2 && j > (Ncc/2)) {
				in[i*Ncc+j] = in[i*Ncc+Ncc+1-j];
			}
			if(i > (Ncc/2) && j <= (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc+1-i)*Ncc+j];
			}

			if(i > (Ncc/2) && j > (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc+1-i)*Ncc+Ncc+1-j];
			}
		}
	}
}
//--------------------------------------------------------------------
void kernel_shears_iso(int Ncc,REAL *in1,REAL *in2,REAL Dcell) {
	int i,j;
	REAL x,y,r;
	REAL epsilon = 0.00000001*Dcell;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <=(Ncc/2)  && j <=(Ncc/2)) {
			x = (REAL)(i)*Dcell-0.5*Dcell;
			y = (REAL)(j)*Dcell-0.5*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);

			if(r > Dcell*(REAL)Ncc/2.0) {
				in1[i*Ncc+j] = 0.0;
				in2[i*Ncc+j] = 0.0;
			}
			else {
				in1[i*Ncc+j] =  (y*y-x*x)/(M_PI*r*r*r*r);
				in2[i*Ncc+j] = (-2.0*x*y)/(M_PI*r*r*r*r);
			}
		}

		else {
			if(i <= (Ncc/2) && j > (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[i*Ncc+Ncc+1-j];
				in2[i*Ncc+j]  = -in2[i*Ncc+Ncc+1-j];
			}
			if(i > (Ncc/2) && j <= (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[(Ncc+1-i)*Ncc+j];
				in2[i*Ncc+j]  = -in2[(Ncc+1-i)*Ncc+j];
			}

			if(i > (Ncc/2) && j > (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[(Ncc+1-i)*Ncc+Ncc+1-j];
				in2[i*Ncc+j]  =  in2[(Ncc+1-i)*Ncc+Ncc+1-j];
			}
		}
	}
}
//--------------------------------------------------------------------
void kernel_alphas_iso(int Ncc,REAL *in1,REAL *in2,REAL Dcell) {
	int i,j;
	REAL x,y,r;
	REAL epsilon = 0.00000001*Dcell;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <=(Ncc/2)  && j <=(Ncc/2)) {
			x = (REAL)(i)*Dcell;//-0.5*Dcell;
			y = (REAL)(j)*Dcell;//-0.5*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);

			if(r > Dcell*(REAL)Ncc/2.0) {
				in1[i*Ncc+j] = 0.0;
				in2[i*Ncc+j] = 0.0;
			}
			else {
				in1[i*Ncc+j] = x/(M_PI*r*r);
				in2[i*Ncc+j] = y/(M_PI*r*r);
			}

		}
		else {
			if(i <= Ncc/2 && j > (Ncc/2)) {
				in1[i*Ncc+j]  =  in1[i*Ncc+Ncc+1-j];
				in2[i*Ncc+j]  = -in2[i*Ncc+Ncc+1-j];
			}
			if(i > (Ncc/2) && j <= (Ncc/2)) {
				in1[i*Ncc+j]  = -in1[(Ncc+1-i)*Ncc+j];
				in2[i*Ncc+j]  =  in2[(Ncc+1-i)*Ncc+j];
			}

			if(i > (Ncc/2) && j > (Ncc/2)) {
				in1[i*Ncc+j]  = -in1[(Ncc+1-i)*Ncc+Ncc+1-j];
				in2[i*Ncc+j]  = -in2[(Ncc+1-i)*Ncc+Ncc+1-j];
			}
		}
	}
}
//--------------------------------------------------------------------
void kernel_smooth_iso(double sigma,int Ncc,REAL *in,REAL Dcell) {
	int i,j;
	REAL x,y,r;
	REAL epsilon = 0.00000001*Dcell;
	REAL cnorm = 0.0;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if(i <=(Ncc/2)  && j <=(Ncc/2)) {
			x = (REAL)(i)*Dcell-0.5*Dcell;
			y = (REAL)(j)*Dcell-0.5*Dcell;
			r = sqrt(x*x+y*y+epsilon*epsilon);

			in[i*Ncc+j] = 1.0/(2.0*M_PI*sigma*sigma)*exp(-(r*r)/(2.0*sigma*sigma));
		}
		else {
			if(i <= Ncc/2 && j > (Ncc/2)) {
				in[i*Ncc+j] = in[i*Ncc+Ncc+1-j];
			}
			if(i > (Ncc/2) && j <= (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc+1-i)*Ncc+j];
			}

			if(i > (Ncc/2) && j > (Ncc/2)) {
				in[i*Ncc+j] = in[(Ncc+1-i)*Ncc+Ncc+1-j];
			}
		}
		cnorm += in[i*Ncc+j]*Dcell*Dcell;
	}

	double ctotal = 0.0;
	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		in[i*Ncc+j] = in[i*Ncc+j]/cnorm;
		ctotal += in[i*Ncc+j]*Dcell*Dcell;
	}
}
//--------------------------------------------------------------------
void lanczos_diff_1_tag(REAL *mi, REAL *m1, REAL *m2, REAL Dcell, int Ncc, int dif_tag) {
	long i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
	long index;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if (i==0) {i_m1 = Ncc-1;i_m2 = Ncc-2;i_m3 = Ncc-3;}
		else if (i==1) {i_m1 = 0;i_m2 = Ncc-1;i_m3 = Ncc-2;}
		else if (i==2) {i_m1 = 1;i_m2 = 0;i_m3 = Ncc-1;}
		else {i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;}
		if (j==0) {j_m1 = Ncc-1;j_m2 = Ncc-2;j_m3 = Ncc-3;}
		else if (j==1) {j_m1 = 0;j_m2 = Ncc-1;j_m3 = Ncc-2;}
		else if (j==2) {j_m1 = 1;j_m2 = 0;j_m3 = Ncc-1;}
		else {j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;}
		if (i==Ncc-1) {i_p1 = 0;i_p2 = 1;i_p3 = 2;}
		else if (i==Ncc-2) {i_p1 = Ncc-1;i_p2 = 0;i_p3 = 1;}
		else if (i==Ncc-3) {i_p1 = Ncc-2;i_p2 = Ncc-1;i_p3 = 0;}
		else {i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;}
		if (j==Ncc-1) {j_p1 = 0;j_p2 = 1;j_p3 = 2;}
		else if (j==Ncc-2) {j_p1 = Ncc-1;j_p2 = 0;j_p3 = 1;}
		else if (j==Ncc-2) {j_p1 = Ncc-2;j_p2 = Ncc-1;j_p3 = 0;}
		else {j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;}

		index = i*Ncc+j;


		if (dif_tag==0) {
			m1[index] = (mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])*2.0/3.0/Dcell
			    	  - (mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])/12.0/Dcell;
			m2[index] = (mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])*2.0/3.0/Dcell
					  - (mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])/12.0/Dcell;
		}

		if (dif_tag==1) {
			m1[index] =(1.0*(mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])
			  		  + 2.0*(mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])
			   		  + 3.0*(mi[i_p3*Ncc+j]-mi[i_m3*Ncc+j]))/(28.0*Dcell);
			m2[index] =(1.0*(mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])
			  		  + 2.0*(mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])
			   		  + 3.0*(mi[i*Ncc+j_p3]-mi[i*Ncc+j_m3]))/(28.0*Dcell);
		}

		if (dif_tag==2) {
			m1[index] =(5.0*(mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])
			  	   	  + 4.0*(mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])
			     	  + 1.0*(mi[i_p3*Ncc+j]-mi[i_m3*Ncc+j]))/(32.0*Dcell);
			m2[index] =(5.0*(mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])
			  		  + 4.0*(mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])
			  		  + 1.0*(mi[i*Ncc+j_p3]-mi[i*Ncc+j_m3]))/(32.0*Dcell);
		}

		if (dif_tag==3) {
			m1[index] = (58.0*(mi[i_p1*Ncc+j]-mi[i_m1*Ncc+j])
			            + 67.0*(mi[i_p2*Ncc+j]-mi[i_m2*Ncc+j])
			  	   	 	+ 22.0*(mi[i_p3*Ncc+j]-mi[i_m3*Ncc+j]))/(252.0*Dcell);
			m2[index] = (58.0*(mi[i*Ncc+j_p1]-mi[i*Ncc+j_m1])
					    + 67.0*(mi[i*Ncc+j_p2]-mi[i*Ncc+j_m2])
						- 22.0*(mi[i*Ncc+j_p3]-mi[i*Ncc+j_m3]))/(252.0*Dcell);
		}
	}
}
//--------------------------------------------------------------------
void lanczos_diff_2_tag(REAL *m1, REAL *m2, REAL *m11, REAL *m12, REAL *m21, REAL *m22, REAL Dcell, int Ncc, int dif_tag) {
	long i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
	long index;

	for(i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		if (i==0) {i_m1 = Ncc-1;i_m2 = Ncc-2;i_m3 = Ncc-3;}
		else if (i==1) {i_m1 = 0;i_m2 = Ncc-1;i_m3 = Ncc-2;}
		else if (i==2) {i_m1 = 1;i_m2 = 0;i_m3 = Ncc-1;}
		else {i_m1 = i-1;i_m2 = i-2;i_m3 = i-3;}
		if (j==0) {j_m1 = Ncc-1;j_m2 = Ncc-2;j_m3 = Ncc-3;}
		else if (j==1) {j_m1 = 0;j_m2 = Ncc-1;j_m3 = Ncc-2;}
		else if (j==2) {j_m1 = 1;j_m2 = 0;j_m3 = Ncc-1;}
		else {j_m1 = j-1;j_m2 = j-2;j_m3 = j-3;}
		if (i==Ncc-1) {i_p1 = 0;i_p2 = 1;i_p3 = 2;}
		else if (i==Ncc-2) {i_p1 = Ncc-1;i_p2 = 0;i_p3 = 1;}
		else if (i==Ncc-3) {i_p1 = Ncc-2;i_p2 = Ncc-1;i_p3 = 0;}
		else {i_p1 = i+1;i_p2 = i+2;i_p3 = i+3;}
		if (j==Ncc-1) {j_p1 = 0;j_p2 = 1;j_p3 = 2;}
		else if (j==Ncc-2) {j_p1 = Ncc-1;j_p2 = 0;j_p3 = 1;}
		else if (j==Ncc-2) {j_p1 = Ncc-2;j_p2 = Ncc-1;j_p3 = 0;}
		else {j_p1 = j+1;j_p2 = j+2;j_p3 = j+3;}

		index = i*Ncc+j;

		if (dif_tag==0) {
			m11[index] = (m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])*2.0/3.0/Dcell
			  		   - (m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])/12.0/Dcell;
			m22[index] = (m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])*2.0/3.0/Dcell
			  		   - (m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])/12.0/Dcell;
			m12[index] = (m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])*2.0/3.0/Dcell
			  		   - (m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])/12.0/Dcell;
			m21[index] = (m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])*2.0/3.0/Dcell
					   - (m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])/12.0/Dcell;
		}

		if (dif_tag==1) {
			m11[index] =(1.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
			  		   + 2.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
			   		   + 3.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(28.0*Dcell);
			m22[index] =(1.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
			  		   + 2.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
			   		   + 3.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(28.0*Dcell);
			m12[index] =(1.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
			  		   + 2.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
			   		   + 3.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(28.0*Dcell);
			m21[index] =(1.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
					   + 2.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
			     	   + 3.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(28.0*Dcell);
		}

		if (dif_tag==2) {
			m11[index] = (5.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
			  	   	    + 4.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
			     	 	+ 1.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(32.0*Dcell);
			m22[index] = (5.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
			  		    + 4.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
			  		    + 1.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(32.0*Dcell);
			m12[index] = (5.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
			  		    + 4.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
			  		    + 1.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(32.0*Dcell);
			m21[index] = (5.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
					    + 4.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
					 	+ 1.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(32.0*Dcell);
		}

		if (dif_tag==3) {
			m11[index] = (58.0*(m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])
			            + 67.0*(m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])
			  	   	 	+ 22.0*(m1[i_p3*Ncc+j]-m1[i_m3*Ncc+j]))/(252.0*Dcell);
			m22[index] = (58.0*(m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])
					    + 67.0*(m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])
					    - 22.0*(m2[i*Ncc+j_p3]-m2[i*Ncc+j_m3]))/(252.0*Dcell);
			m12[index] = (58.0*(m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])
					    + 67.0*(m1[i*Ncc+j_p2]-m1[i*Ncc+j_m2])
						- 22.0*(m1[i*Ncc+j_p3]-m1[i*Ncc+j_m3]))/(252.0*Dcell);
			m21[index] = (58.0*(m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])
					    + 67.0*(m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])
					 	- 22.0*(m2[i_p3*Ncc+j]-m2[i_m3*Ncc+j]))/(252.0*Dcell);
		}
	}
}
//--------------------------------------------------------------------
void write_2_signals(char *out1,char *out2,REAL *in1,REAL *in2, int Ncc) {
	long i,j;
	long index;
	FILE *f1,*f2;

	f1 =fopen(out1,"wb");
	f2 =fopen(out2,"wb");

	for (i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		index = i*Ncc+j;

		fwrite(&in1[index],sizeof(REAL),1,f1);
		fwrite(&in2[index],sizeof(REAL),1,f2);
	}
	fclose(f1);
	fclose(f2);
}
//--------------------------------------------------------------------
void write_1_signal(char *out1,REAL *in1, int Ncc) {
	long i,j;
	long index;
	FILE *f1;

	f1 =fopen(out1,"wb");

	for (i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		fwrite(&in1[index],sizeof(REAL),1,f1);
	}
	fclose(f1);
}
//--------------------------------------------------------------------
void print_1_signal(REAL *in1,int Ncc) {
	long i,j;
	long index;

	for (i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		printf("%lf ",in1[index]);
	}
}
//--------------------------------------------------------------------
void print_2_signals(REAL *in1,REAL *in2,int Ncc) {
	long i,j;
	long index;

	for (i=0;i<Ncc;i++) for(j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		printf("%lf %lf ",in1[index],in2[index]);
	}
}
//--------------------------------------------------------------------
void kappa_to_phi(REAL *kappa, REAL *phi,int Ncc2,REAL Dcell) {

	REAL *green_iso = (REAL *)calloc(Ncc2*Ncc2,sizeof(REAL));
	kernel_green_iso(Ncc2,green_iso,Dcell);
	REAL *phi_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));
	convolve_fft(kappa,green_iso,phi_tmp,Ncc2,Ncc2,Dcell,Dcell);
	corner_matrix(phi_tmp,Ncc2,Ncc2,phi);
	free(phi_tmp);
}
//--------------------------------------------------------------------
void kappa_to_alphas(REAL *kappa,REAL *alpha1,REAL *alpha2,int Ncc2,REAL Dcell) {

	REAL *alpha1_iso = (REAL *)calloc(Ncc2*Ncc2,sizeof(REAL));
	REAL *alpha2_iso = (REAL *)calloc(Ncc2*Ncc2,sizeof(REAL));

	kernel_alphas_iso(Ncc2,alpha1_iso,alpha2_iso,Dcell);

	REAL *alpha1_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));
	REAL *alpha2_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));

	convolve_fft(kappa,alpha1_iso,alpha1_tmp,Ncc2,Ncc2,Dcell,Dcell);
	convolve_fft(kappa,alpha2_iso,alpha2_tmp,Ncc2,Ncc2,Dcell,Dcell);
	free(alpha1_iso);
	free(alpha2_iso);

	corner_matrix(alpha1_tmp,Ncc2,Ncc2,alpha1);
	corner_matrix(alpha2_tmp,Ncc2,Ncc2,alpha2);

	free(alpha1_tmp);
	free(alpha2_tmp);
}
//--------------------------------------------------------------------
void kappa_to_shears(REAL *kappa,REAL *shear1,REAL *shear2, int Ncc2,REAL Dcell) {

	REAL *shear1_iso = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));
	REAL *shear2_iso = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));

	kernel_shears_iso(Ncc2,shear1_iso,shear2_iso,Dcell);

	REAL *shear1_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));
	REAL *shear2_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));

	convolve_fft(kappa,shear1_iso,shear1_tmp,Ncc2,Ncc2,Dcell,Dcell);
	convolve_fft(kappa,shear2_iso,shear2_tmp,Ncc2,Ncc2,Dcell,Dcell);

	corner_matrix(shear1_tmp,Ncc2,Ncc2,shear1);
	corner_matrix(shear2_tmp,Ncc2,Ncc2,shear2);

	free(shear1_tmp);
	free(shear2_tmp);
	free(shear1_iso);
	free(shear2_iso);
}
//--------------------------------------------------------------------
void kappa_to_kappac(REAL *kappa,REAL *kappa_cal, int Ncc2,REAL Dcell) {

	REAL *shear1_iso = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));
	REAL *shear2_iso = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));

	kernel_shears_iso(Ncc2,shear1_iso,shear2_iso,Dcell);

	REAL *shear1_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));
	REAL *shear2_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));

	convolve_fft(kappa,shear1_iso,shear1_tmp,Ncc2,Ncc2,Dcell,Dcell);
	convolve_fft(kappa,shear2_iso,shear2_tmp,Ncc2,Ncc2,Dcell,Dcell);

	REAL *kappa1_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));
	REAL *kappa2_tmp = (REAL *)malloc(Ncc2*Ncc2*sizeof(REAL));

	convolve_fft(shear1_tmp,shear1_iso,kappa1_tmp,Ncc2,Ncc2,Dcell,Dcell);
	convolve_fft(shear2_tmp,shear2_iso,kappa2_tmp,Ncc2,Ncc2,Dcell,Dcell);

	corner_matrix(kappa1_tmp,Ncc2,Ncc2,kappa_cal);

	free(shear1_tmp);
	free(shear2_tmp);
	free(shear1_iso);
	free(shear2_iso);
	free(kappa1_tmp);
	free(kappa2_tmp);
}
//--------------------------------------------------------------------
void sdens_to_phi(REAL * sdens, int Nc, REAL boxsize, REAL zl, REAL zs, REAL * phi) {

	int Nc2 = Nc*2;
	REAL dsx = boxsize/(REAL)Nc;

	REAL *kappa0 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	sdens_to_kappa(sdens, Nc, dsx, zl, zs, kappa0);
	REAL *kappa = (REAL *)calloc(Nc2*Nc2,sizeof(REAL));
	zero_padding(kappa0,Nc,Nc,kappa);
	free(kappa0);
	kappa_to_phi(kappa, phi, Nc2, dsx);
	free(kappa);
}
//--------------------------------------------------------------------
void sdens_to_kappac(REAL * sdens, int Nc, REAL boxsize, REAL zl, REAL zs, REAL * kappac) {

	int Nc2 = Nc*2;
	REAL dsx = boxsize/(REAL)Nc;

	REAL *kappa0 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	sdens_to_kappa(sdens, Nc, dsx, zl, zs, kappa0);
	REAL *kappa = (REAL *)calloc(Nc2*Nc2,sizeof(REAL));
	zero_padding(kappa0,Nc,Nc,kappa);
	free(kappa0);
	kappa_to_kappac(kappa, kappac, Nc2, dsx);
	free(kappa);
}
//--------------------------------------------------------------------
void sdens_to_alphas(REAL * sdens, int Nc, REAL boxsize, REAL zl, REAL zs, REAL * alpha1, REAL * alpha2) {

	int Nc2 = Nc*2;
	REAL dsx = boxsize/(REAL)Nc;
	REAL *kappa0 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	sdens_to_kappa(sdens, Nc, dsx, zl, zs, kappa0);
	REAL *kappa = (REAL *)calloc(Nc2*Nc2,sizeof(REAL));
	zero_padding(kappa0,Nc,Nc,kappa);
	free(kappa0);
	kappa_to_alphas(kappa,alpha1,alpha2,Nc2,dsx);
	free(kappa);
}
//--------------------------------------------------------------------
void nkappa_to_alphas(REAL * nkappa, int Nc, REAL boxsize, REAL * alpha1, REAL * alpha2) {

	int Nc2 = Nc*2;
	REAL dsx = boxsize/(REAL)Nc;

	REAL *kappa = (REAL *)calloc(Nc2*Nc2,sizeof(REAL));
	zero_padding(nkappa,Nc,Nc,kappa);
	kappa_to_alphas(kappa,alpha1,alpha2,Nc2,dsx);
	free(kappa);
}
//--------------------------------------------------------------------
void sdens_to_shears(REAL * sdens, int Nc, REAL boxsize, REAL zl, REAL zs, REAL * shear1, REAL * shear2) {


	int Nc2 = Nc*2;
	REAL dsx = boxsize/(REAL)Nc;

	REAL *kappa0 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	sdens_to_kappa(sdens, Nc, dsx, zl, zs, kappa0);
	REAL *kappa = (REAL *)calloc(Nc2*Nc2,sizeof(REAL));
	zero_padding(kappa0,Nc,Nc,kappa);
	free(kappa0);
	kappa_to_shears(kappa,shear1,shear2,Nc2,dsx);
	free(kappa);

	//REAL *kappa = calloc(Nc2*Nc2,sizeof(REAL));
	//sdens_to_kappa(pmass, sdens, Nc, dsx, zl, zs, kappa);
	//kappa_to_shears(kappa,shear1,shear2,Nc2,dsx);
	//free(kappa);
}
//--------------------------------------------------------------------
void sdens_to_mu(REAL *sdens, int Nc, REAL boxsize, REAL zl, REAL zs, REAL *mu) {

	REAL *kappac = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	sdens_to_kappac(sdens, Nc, boxsize, zl, zs, kappac);

	REAL *shear1 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	REAL *shear2 = (REAL *)calloc(Nc*Nc,sizeof(REAL));
	sdens_to_shears(sdens, Nc, boxsize, zl, zs, shear1, shear2);

	int i,j,index;
	for (i=0;i<Nc;i++) for (j=0;j<Nc;j++) {
		index = i*Nc+j;
		mu[index] = 1.0/((1.0-kappac[index])*(1.0-kappac[index])
				-shear1[index]*shear1[index]-shear2[index]*shear2[index]);
	}

	free(kappac);
	free(shear1);
	free(shear2);
}
//--------------------------------------------------------------------
void alphas_to_mu(REAL *alpha1,REAL *alpha2,int Ncc, REAL Dcell, REAL *mu) {

	int i,j,index;
	REAL *al11 = (REAL *)malloc(Ncc*Ncc*sizeof(REAL));
	REAL *al12 = (REAL *)malloc(Ncc*Ncc*sizeof(REAL));
	REAL *al21 = (REAL *)malloc(Ncc*Ncc*sizeof(REAL));
	REAL *al22 = (REAL *)malloc(Ncc*Ncc*sizeof(REAL));

	lanczos_diff_2_tag(alpha1,alpha2,al11, al12, al21, al22, Dcell, Ncc,1);
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		mu[index] = 1.0/(1.0-(al11[index]+al22[index])+al11[index]*al22[index]-al12[index]*al21[index]);
	}
	free(al11);
	free(al12);
	free(al21);
	free(al22);
}
//--------------------------------------------------------------------
void calculate_mu(REAL *kappa,REAL *shear1,REAL *shear2, int Ncc, REAL *mu) {

	int i,j,index;
	for (i=0;i<Ncc;i++) for (j=0;j<Ncc;j++) {
		index = i*Ncc+j;
		mu[index] = 1.0/((1.0-kappa[index])*(1.0-kappa[index])
				-shear1[index]*shear1[index]-shear2[index]*shear2[index]);
	}
}
//--------------------------------------------------------------------
void find_critical_curves(REAL *mu_in, int nx, int ny,REAL *mu_out) {

    REAL mu_ip1 = 0.0;
    REAL mu_im1 = 0.0;
    REAL mu_jp1 = 0.0;
    REAL mu_jm1 = 0.0;
    REAL sign_total = 0.0;

	int im1, jm1, ip1, jp1;
	int i, j, index;
    for (i = 0;i < nx;i++) for (j = 0;j < ny;j++) {
		index = i*ny+j;

		if (i==0) im1 = nx-1;
		else im1 = i-1;
		if (j==0) jm1 = ny-1;
		else jm1 = j-1;
		if (i==nx-1) ip1 = 0;
		else ip1 = i+1;
		if (j==ny-1) jp1 = 0;
		else jp1 = j+1;

		mu_ip1 = mu_in[ip1*ny+j];
		mu_im1 = mu_in[im1*ny+j];
		mu_jp1 = mu_in[i*ny+jp1];
		mu_jm1 = mu_in[i*ny+jm1];

		sign_total = sign(mu_in[index])
		    *(sign(mu_ip1)+sign(mu_im1)
		     +sign(mu_jp1)+sign(mu_jm1));


		//if (sign_total < 4.0) {
		//	rcc_tmp = sqrt(posx1[index]*posx1[index] + posx2[index]*posx2[index]);
		//	rcc_max = (rcc_tmp>=rcc_max)? rcc_tmp:rcc_max;
		//}
		if (sign_total < 4.0) {
			mu_out[index] = 1.0;
		}
	}
}
