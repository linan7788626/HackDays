__kernel void nie_alphas_cl(                      
   __const __global float* input1,                  
   __const __global float* input2,                  
   __const __global float* lpar,                  
   __global float* output1,                 
   __global float* output2,                 
   const int count) {                                          
	int i = get_global_id(0);               


    local float xc1;
	xc1 = lpar[0];
    local float xc2;
	xc2 = lpar[1];
    local float q;
	q   = lpar[2];
    local float rc;
	rc  = lpar[3];
	local float re;
    re  = lpar[4];
	local float pha;
    pha = lpar[5];

	local float xi1;
	xi1 = input1[i];
	local float xi2;
	xi2 = input2[i];

    local float phirad;
	phirad = pha*3.141592653589793/180.0;;
    local float cosa;
	cosa = cos(phirad);
    local float sina; 
	sina = sin(phirad);
	local float phi,a1,a2;
	local float xt1,xt2;

	xt1 = (xi1-xc1)*cosa+(xi2-xc2)*sina;
	xt2 = (xi2-xc2)*cosa-(xi1-xc1)*sina;
	phi = sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc);

	a1 = sqrt(q)/sqrt(1.0-q*q)*atan(sqrt(1.0-q*q)*xt1/(phi+rc/q));
	a2 = sqrt(q)/sqrt(1.0-q*q)*atanh(sqrt(1.0-q*q)*xt2/(phi+rc*q));

	output1[i] = (a1*cosa-a2*sina)*re;
	output2[i] = (a2*cosa+a1*sina)*re;
}
