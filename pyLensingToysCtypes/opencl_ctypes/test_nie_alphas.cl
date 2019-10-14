__kernel void nie_alphas_cl(                      
   __const __global float2 * inputs,                  
   __const __global float * lpar,                  
   __global float2 * outputs,                 
   const int count) {                                          

	int i = get_global_id(0);               

    private float xc1;
	xc1 = lpar[0];
    private float xc2;
	xc2 = lpar[1];
    private float q;
	q   = lpar[2];
    private float rc;
	rc  = lpar[3];
	private float re;
    re  = lpar[4];
	private float pha;
    pha = lpar[5];

	private float xi1;
	xi1 = inputs[i].x;
	private float xi2;
	xi2 = inputs[i].y;

    private float phirad;
	phirad = pha*3.141592653589793f/180.0f;;
    private float cosa;
	cosa = cos(phirad);
    private float sina; 
	sina = sin(phirad);
	private float phi,a1,a2;
	private float xt1,xt2;

	xt1 = (xi1-xc1)*cosa+(xi2-xc2)*sina;
	xt2 = (xi2-xc2)*cosa-(xi1-xc1)*sina;
	phi = sqrt(xt2*xt2+xt1*q*xt1*q+rc*rc);

	a1 = sqrt(q)/sqrt(1.0f-q*q)*atan(sqrt(1.0f-q*q)*xt1/(phi+rc/q));
	a2 = sqrt(q)/sqrt(1.0f-q*q)*atanh(sqrt(1.0f-q*q)*xt2/(phi+rc*q));

	outputs[i].x = (a1*cosa-a2*sina)*re;
	outputs[i].y = (a2*cosa+a1*sina)*re;
}
