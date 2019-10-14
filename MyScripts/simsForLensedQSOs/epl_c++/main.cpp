/* complex.cpp */
/*
  Program to test out C++ complex class
*/

#include <complex.h>
#include <stdio.h>

/* Define complex double data type */
typedef complex<double> dcomp;

int main()
{
  dcomp i, a, b, c, d, e, p, q, r; // Declare complex double variables
  double x, y;

  /* Set complex double variable equal to complex double constant */
  i = dcomp (0., 1.);
  printf("\ni = (%6.4f, %6.4f)\n", i);

  /* Test arithmetic operations with complex double variables */
  a = i * i;
  b = 1. / i;
  printf("\ni*i = (%6.4f, %6.4f)\n", a);
  printf("1/i = (%6.4f, %6.4f)\n", b);

  /* Test mathematical functions using complex double variables */
  c = sqrt(i);
  d = sin(i);
  e = pow(i, 0.25);
  printf("\nsqrt(i) = (%6.4f, %6.4f)\n", c);
  printf("sin(i) = (%6.4f, %6.4f)\n", d);
  printf("i^0.25 = (%6.4f, %6.4f)\n", e);

  /* Test complex operations */
  p = conj(i);
  q = real(i);
  r = imag(i);
  printf("\nconj(i) = (%6.4f, %6.4f)\n", p);
  printf("real(i) = %6.4f\n", q);
  printf("imag(i) = %6.4f\n", r);

  return 0;
}
