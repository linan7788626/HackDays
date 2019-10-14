// Illustration of compiler <complex> class.
#include <iostream>
#include <complex>
using namespace std;

typedef complex<double> dcmplx;

int main(){

  dcmplx a(5.,6.),b;

  cout << "Enter b: ";
  cin >> b;

  cout << "a = " << a << "\n";
  cout << "b = " << b << "\n";

  cout << "a + b = " << a + b << "\n";
  cout << "a * b = " << a * b << "\n";
  cout << "a / b = " << a / b << "\n";
  cout << "|a| = "   << abs(a) << "\n";
  cout << "complex conjugate of a = " << conj(a) << "\n";
  cout << "norm of a = " << norm(a) << "\n";
  cout << "abs of a = " << abs(a) << "\n";
  cout << "exp(a) = " << exp(a) << "\n";
}
