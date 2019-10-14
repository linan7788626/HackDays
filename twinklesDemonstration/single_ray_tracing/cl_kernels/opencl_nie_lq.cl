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
