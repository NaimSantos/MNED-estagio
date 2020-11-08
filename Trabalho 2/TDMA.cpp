#include <iostream>
#include <vector>
#include "alg_utilities.h"

int main(int argc, char **argv) {
	
	size_t N = 11;
	double delta_x = static_cast<double>(1.0/N);
	double delta_t = 0.001;
	double r = delta_t/(delta_x*delta_x);

	std::vector<double> a(N-1, -r/2.0);
	std::vector<double> b(N, 1.0+r);
	std::vector<double> c(N-1, -r/2.0);
	std::vector<double> d(N, 0.0);
	std::vector<double> f(N, 0.0);

	f[5] = 1; f[6] = 2; f[7] = 1;

	std::cout << "f = (";
	for (size_t i=0; i<N; i++) {
		std::cout << f[i];
		if (i < N-1) {
			std:: cout << ", ";
		}
	}
	std::cout << ")" << std::endl << std::endl;
	for (size_t i=1; i<N-1; i++) {
		d[i] = r*0.5*f[i+1] + (1.0-r)*f[i] + r*0.5*f[i-1];
	}

	tdma_solver(a, b, c, d, f);

	std::cout << "f = (";
	for (size_t i=0; i<N; i++) {
		std::cout << f[i];
		if (i < N-1) {
			std:: cout << ", ";
		}
	}
	std::cout << ")" << std::endl;

	return 0;
}
