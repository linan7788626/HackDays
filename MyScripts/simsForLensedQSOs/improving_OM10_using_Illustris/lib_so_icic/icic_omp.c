#include <stdio.h>
#include <stdlib.h>

void inverse_cic(double *source_map, double *posy1, double *posy2, double ysc1, double ysc2,double dsi, int nsx, int nsy, int nlimgs, double *lensed_map) {

	int i1,j1,i;
	double xb1,xb2;
	double ww1,ww2,ww3,ww4,wx,wy;

	for(i=0;i<nlimgs;i++){

		xb1 = (posy1[i]-ysc1)/dsi+(double)nsx/2.0-0.5;
		xb2 = (posy2[i]-ysc2)/dsi+(double)nsy/2.0-0.5;

		i1 = (int)xb1;
		j1 = (int)xb2;

		wx = 1.-(xb1-(double)(i1));
		wy = 1.-(xb2-(double)(j1));

		ww1 = wx*wy;
		ww2 = wx*(1.0-wy);
		ww3 = (1.0-wx)*wy;
		ww4 = (1.0-wx)*(1.0-wy);

		if (i1<0||i1>nsx-2||j1<0||j1>nsy-2) continue;

		lensed_map[i] = ww1*source_map[i1*nsx+j1]
						  + ww2*source_map[i1*nsx+j1+1]
						  + ww3*source_map[(i1+1)*nsx+j1]
						  + ww4*source_map[(i1+1)*nsx+j1+1];
	}
}
