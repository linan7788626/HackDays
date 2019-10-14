#include <iostream>
#include <complex>
#include <gsl/gsl_sf_hyperg.h>

typedef std::complex<double> dcmplx;


#include <boost/numeric/odeint.hpp>

namespace {

	std::array<dcmplx, 2> hypser(dcmplx const &A, dcmplx const &B, dcmplx const &C, dcmplx const &Z) {
		dcmplx series, deriv, fac{1};
		auto temp{fac}, a{A}, b{B}, c{C};
		for (std::size_t n{1}; n <= 1000; ++n) {
			fac *= a * b / c;
			deriv += fac;
			fac *= 1. / n * Z;
			series = temp + fac;
			if (series == temp)
				return {series, deriv};
			temp = series;
			a += 1;
			b += 1;
			c += 1;
		}
		throw("convergence failure in hypser");
	}

} // namespace

dcmplx hyp2f1(dcmplx const &A, dcmplx const &B, dcmplx const &C, dcmplx const &Z) {
	dcmplx z0;
	if (std::norm(Z) <= 0.25)
		return hypser(A, B, C, Z).front();
	else if (std::real(Z) < 0)
		z0 = -0.5;
	else if (std::real(Z) <= 1)
		z0 = 0.5;
	else {
		using namespace std::literals;
		z0 = std::imag(Z) >= 0 ? 0.5i : -0.5i;
	}
	auto const Dz{Z - z0};
	auto dependentVariable{hypser(A, B, C, z0)};
	boost::numeric::odeint::integrate_adaptive(
	    boost::numeric::odeint::bulirsch_stoer<decltype(dependentVariable)>{1e-14, 1e-14},
	    [&](auto const &DependentVariable, auto &derivative, double const Variable)
	    {
			derivative.front() = DependentVariable.back() * Dz;
			auto const z{z0 + Variable * Dz};
			derivative.back() = (A * B * DependentVariable.front() - (C -
						(A + B + 1.) * z) * DependentVariable.back()) * Dz
						/ (z * (1. - z));
	    },
	    dependentVariable, 0., 1., 0.1);
	return dependentVariable.front();
}



double kappa_epl(double x1, double x2, double b, double q, double t) {
	double R, res;
	R = sqrt(q*q*x1*x1+x2*x2);
	res = 0.5*(2.0-t)*pow((b/R),t);
	return res;
}

dcmplx alphas_epl(double x1, double x2, double b, double q, double t) {
	dcmplx alphas;
	//polar (2.0, 0.5);
	double R, phi;
	R = sqrt(q*q*x1*x1+x2*x2);
	phi = atan2(q*x1,x2);
	dcmplx e_phi;
	e_phi = std::polar(1.0, phi);
	dcmplx e_2phi;
	e_2phi = std::polar(1.0, phi*2.0);

	dcmplx res;


	double fa = 1.0;
	double fb = 0.5*t;
	double fc = 2.0 - 0.5*t;
	dcmplx fd = -(1.0-q)/(1.0+q)*e_2phi;

	dcmplx F_2f1 = hyp2f1(fa, fb, fc, fd);

	alphas = 2.0*b/(1.0+q)*pow((b/R), t-1)*e_phi*F_2f1;
	return alphas;
}

int main() {

	dcmplx a(5.0,6.0);
	a = alphas_epl(0.5, 0.6, 2.0, 0.7, 1.2);
	dcmplx b(0.1,0.4);

	std::cout << "a = " << a << "\n";
	std::cout << "b = " << b << "\n";
}
