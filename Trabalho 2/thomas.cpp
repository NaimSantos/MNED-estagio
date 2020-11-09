#include <iostream>
#include <vector>
#include <fstream>
#include "alg_utilities.h"

int main(int argc, char **argv) {

	//MFD para equação do calor. 15 nós
	const unsigned int N = 15;
	double delta_x = 1.0/N;
	double delta_t = 0.001;
	double r = delta_t/(delta_x*delta_x);

	//Creação de vetores para armazenar os coeficientes
	std::vector<double> a(N-1, -r/2.0);
	std::vector<double> b(N, 1.0+r);
	std::vector<double> c(N-1, -r/2.0);
	std::vector<double> d(N, 0.0);
	std::vector<double> f(N, 0.0);

	//Preencher o valor inicial do vetor passo de tempo atual usando 3 picos 
	f[5] = 1; f[6] = 2; f[7] = 1;

	// Exibir o vetor solução f antes de cada novo passo de tempo ser calculado
	std::cout << "f = (";
	for (int i=0; i<N; i++) {
		std::cout << f[i];
		if (i < N-1) {
			std:: cout << ", ";
		}
	}
	std::cout << ")" << std::endl << std::endl;

	//Preencher o vetor d atual (tempos)
	for (int i=1; i<N-1; i++) {
		d[i] = r*0.5*f[i+1] + (1.0-r)*f[i] + r*0.5*f[i-1];
	}

	//Use o TDMA para resolver o sistema tridiagonal:
	thomas_solver(a, b, c, d, f);

	// Exibir o vetor solução:
	std::cout << "f = (";
	for (int i=0; i<N; i++) {
		std::cout << f[i];
		if (i < N-1) {
			std:: cout << ", ";
		}
	}
	std::cout << ")" << std::endl;

	
	return 0;
}
