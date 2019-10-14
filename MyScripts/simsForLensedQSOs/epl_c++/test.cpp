// g++ hg.cc -o hg -L/usr/local/lib -lgsl -lgslcblas -lm && ./hg
#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

// ----------------------------------------------------
gsl_complex operator/(   gsl_complex  z1 , const gsl_complex & z2) {
    return gsl_complex_div(z1, z2);
}
gsl_complex operator*(   gsl_complex  z1 , const gsl_complex & z2) {
    return gsl_complex_mul(z1, z2);
}
gsl_complex operator+(   gsl_complex  z1 , const gsl_complex & z2) {
    return gsl_complex_add(z1, z2);
}
gsl_complex operator+(   gsl_complex z2  , const double  & x ) {
    return gsl_complex_add( z2,  gsl_complex_rect(x,0));
}
gsl_complex operator/(   gsl_complex z2  , const double  & x ) {
    return gsl_complex_div( z2,  gsl_complex_rect(x,0));
}
// ----------------------------------------------------

gsl_complex hy_ge( const gsl_complex & a
,                  const gsl_complex & b
,                  const gsl_complex & c
,                  const gsl_complex & z
,           unsigned int n = 100                   // accuracy
,           unsigned int i = 0                     // itteration step
,            gsl_complex t = gsl_complex_rect(1,0) // coefficient t for internal calculations
) {

    gsl_complex  t_next = (a+i)*(b+i)/(c+i)/(i+1) * t;

    return (i==n+1) ? gsl_complex_rect(0,0) : t*gsl_complex_pow_real(z,i)+hy_ge(a,b,c,z,n,i+1,t_next);
}

// ----------------------------------------------------
int main (void)
{

    gsl_complex a = gsl_complex_rect(  1,  0);
    gsl_complex b = gsl_complex_rect(  2,  0);
    gsl_complex c = gsl_complex_rect(  2,  2);
    gsl_complex z = gsl_complex_rect(0.2,0.1);

    gsl_complex z_out = hy_ge(a,b,c,z);

    printf( " %.10f%+.10f i \n",  GSL_REAL( z_out ), GSL_IMAG( z_out ) ) ;

    return 0;
}
// ----------------------------------------------------
