inline int sign(float x) {
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

inline float deg2rad(float pha) {
	float res = 0;
	res = pha*3.141592653589793/180.0;
	return res;
}

__kernel void lanczos_diff_1_tag(
	__global float *m, 
	__global float *m1, 
	__global float *m2, 
	__global float Dcell, 
	const int Ncc, 
	const int dif_tag) {

    int i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
    int index;

	int i = get_global_id(0);               
	int j = get_global_id(1);               

    if (i<Ncc && j<Ncc) {

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
        if (dif_tag==-1) {
            m1[index] = (m[i_p1*Ncc+j]-m[i_m1*Ncc+j])/(2.0*Dcell);
            m2[index] = (m[i*Ncc+j_p1]-m[i*Ncc+j_m1])/(2.0*Dcell);
        }

        if (dif_tag==0) {
            m1[index] = (m[i_p1*Ncc+j]-m[i_m1*Ncc+j])*2.0/3.0/Dcell
					  - (m[i_p2*Ncc+j]-m[i_m2*Ncc+j])/12.0/Dcell;
            m2[index] = (m[i*Ncc+j_p1]-m[i*Ncc+j_m1])*2.0/3.0/Dcell
					  - (m[i*Ncc+j_p2]-m[i*Ncc+j_m2])/12.0/Dcell;
        }

        if (dif_tag==1) {
            m1[index] =(1.0*(m[i_p1*Ncc+j]-m[i_m1*Ncc+j])
                      + 2.0*(m[i_p2*Ncc+j]-m[i_m2*Ncc+j])
                      + 3.0*(m[i_p3*Ncc+j]-m[i_m3*Ncc+j]))/(28.0*Dcell);
            m2[index] =(1.0*(m[i*Ncc+j_p1]-m[i*Ncc+j_m1])
                      + 2.0*(m[i*Ncc+j_p2]-m[i*Ncc+j_m2])
                      + 3.0*(m[i*Ncc+j_p3]-m[i*Ncc+j_m3]))/(28.0*Dcell);
        }

        if (dif_tag==2) {
            m1[index] = (5.0*(m[i_p1*Ncc+j]-m[i_m1*Ncc+j])
					   + 4.0*(m[i_p2*Ncc+j]-m[i_m2*Ncc+j])
					   + 1.0*(m[i_p3*Ncc+j]-m[i_m3*Ncc+j]))/(32.0*Dcell);
            m2[index] = (5.0*(m[i*Ncc+j_p1]-m[i*Ncc+j_m1])
                       + 4.0*(m[i*Ncc+j_p2]-m[i*Ncc+j_m2])
                       + 1.0*(m[i*Ncc+j_p3]-m[i*Ncc+j_m3]))/(32.0*Dcell);
        }

        if (dif_tag==3) {
            m1[index] = (58.0*(m[i_p1*Ncc+j]-m[i_m1*Ncc+j])
                       + 67.0*(m[i_p2*Ncc+j]-m[i_m2*Ncc+j])
                       + 22.0*(m[i_p3*Ncc+j]-m[i_m3*Ncc+j]))/(252.0*Dcell);
            m2[index] = (58.0*(m[i*Ncc+j_p1]-m[i*Ncc+j_m1])
                       + 67.0*(m[i*Ncc+j_p2]-m[i*Ncc+j_m2])
                       - 22.0*(m[i*Ncc+j_p3]-m[i*Ncc+j_m3]))/(252.0*Dcell);
        }
    }
}

__kernel void lanczos_diff_2_tag(
	__global float *m1, 
	__global float *m2, 
	__global float *m11, 
	__global float *m12, 
	__global float *m21, 
	__global float *m22, 
	__global float Dcell, 
	const int Ncc, 
	const int dif_tag) {

    int i_m3,i_p3,j_m3,j_p3,i_m2,i_p2,j_m2,j_p2,i_m1,j_m1,i_p1,j_p1,i,j;
    int index;

	int i = get_global_id(0);               
	int j = get_global_id(1);               

    if (i<Ncc && j<Ncc) {

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
        if (dif_tag==-1) {
            m11[index] = (m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])/(2.0*Dcell);
            m12[index] = (m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])/(2.0*Dcell);
            m21[index] = (m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])/(2.0*Dcell);
            m22[index] = (m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])/(2.0*Dcell);
        }

        if (dif_tag==0) {
            m11[index] = (m1[i_p1*Ncc+j]-m1[i_m1*Ncc+j])*2.0/3.0/Dcell
            - (m1[i_p2*Ncc+j]-m1[i_m2*Ncc+j])/12.0/Dcell;
            m22[index] = (m2[i*Ncc+j_p1]-m2[i*Ncc+j_m1])*2.0/3.0/Dcell
            - (m2[i*Ncc+j_p2]-m2[i*Ncc+j_m2])/12.0/Dcell;
            m21[index] = (m2[i_p1*Ncc+j]-m2[i_m1*Ncc+j])*2.0/3.0/Dcell
            - (m2[i_p2*Ncc+j]-m2[i_m2*Ncc+j])/12.0/Dcell;
            m12[index] = (m1[i*Ncc+j_p1]-m1[i*Ncc+j_m1])*2.0/3.0/Dcell
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

__kernel void opencl_nie_lq (                      
	__global float* input1,                  
	__global float* input2,                  
	__global float* lpar,                  
	__global float* output1,                 
	__global float* output2,                 
	const int count)               
{                                          
	int i = get_global_id(0);               
	if(i < count) {
		float phirad = pha*3.141592653589793/180.0;
    	float cosa = cos(phirad);
    	float sina = sin(phirad);
		float phi,a1,a2;
		float xt1,xt2;

		xt1 = (input1[i]-xc1)*cosa+(input2[i]-xc2)*sina;
		xt2 = (input2[i]-xc2)*cosa-(input1[i]-xc1)*sina;
		phi = sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc);

		a1 = sqrt(q)/sqrt(1.0-q*q)*atan(sqrt(1.0-q*q)*xt1/(phi+rc/q));
		a2 = sqrt(q)/sqrt(1.0-q*q)*atanh(sqrt(1.0-q*q)*xt2/(phi+rc*q));

		output1[i] = input1[i]-(a1*cosa-a2*sina)*re;
		output2[i] = input2[i]-(a2*cosa+a1*sina)*re;
	}
}

__kernel void opencl_wcic(
	__global float *cic_in,
	__global float *x_in,
	__global float *y_in,
	const float bsx,
	const float bsy,
	const int nx,
	const int ny,
	const int np,
	__global float *cic_out) {

    float dx = bsx/nx;
    float dy = bsy/ny;
    float xc = bsx/2.0;
    float yc = bsy/2.0;
    float wx,wy;
    float xp,yp,zp;

    int i;
    int ip,jp;

	int i = get_global_id(0);               
	if(i < np) {
        xp = (x_in[i]+xc)/dx-0.5;
        yp = (y_in[i]+yc)/dy-0.5;
		zp = cic_in[i];

        ip = (int)xp;
        jp = (int)yp;

		if (ip<0||ip>(nx-2)||jp<0||jp>(ny-2)) break;

        wx = 1.0-(xp-(float)ip);
        wy = 1.0-(yp-(float)jp);

        cic_out[ip*ny+jp] += wx*wy*zp;
        cic_out[ip*ny+(jp+1)] += wx*(1.0-wy)*zp;
        cic_out[(ip+1)*ny+jp] += (1.0-wx)*wy*zp;
        cic_out[(ip+1)*ny+(jp+1)] += (1.0-wx)*(1.0-wy)*zp;
    }
}

__kernel void opencl_gauss_2d(
	__global float *x1,
	__global float *x2,
	__global float *par,
	__global float *res,
	const int count) {

	int i = get_global_id(0);               
	if(i < count) {
		float x1new,x2new;
		float cosa = cos(deg2rad(par[5]));
		float sina = sin(deg2rad(par[5]));

		x1new = (x1[i] - par[0])*cosa+(x2[i]-par[1])*sina);
		x2new = (x2[i] - par[1])*cosa-(x1[i]-par[0])*sina);

		float re_eff = (x1new*x1new)*par[2]+(x2new*x2new)/par[2];
		if (re_eff>4.0) {
			res[i] = 0.0;
		}
		else {
			res[i] = par[3]*exp(-0.5*(re_eff)/(par[4]*par[4]));
		}
	}
}

__kernel void opencl_tophat_2d(
	__global float *x1,
	__global float *x2,
	__global float *par,
	__global float *res,
	const int count) {

	int i = get_global_id(0);               
	if(i < count) {
		float x1new,x2new;
		float cosa = cos(deg2rad(par[5]));
		float sina = sin(deg2rad(par[5]));

		x1new = (x1[i] - par[0])*cosa+(x2[i]-par[1])*sina);
		x2new = (x2[i] - par[1])*cosa-(x1[i]-par[0])*sina);

		float r_ell = sqrt((x1new*x1new)*par[2]+(x2new*x2new)/par[2]);
		if (r_ell>=par[4]) {
			res[i] = -1.0;
		}
		else {
			res[i] = 10000.0;
		}
	}
}

__kernel void opencl_rotate_xy(
	__global float *x1,
	__global float *x2,
	const float xc1,
	const float xc2,
	const float pha,
	__global float *x1_out,
	__global float *x2_out,
	const int count) {

	float cosa = cos(deg2rad(pha));
	float sina = sin(deg2rad(pha));

	int i = get_global_id(0);               
	if(i < count) {
		x1_out[i] = (x1[i] - xc1)*cosa+(x2[i]-xc2)*sina);
		x2_out[i] = (x2[i] - xc2)*cosa-(x1[i]-xc1)*sina);
	}
}
