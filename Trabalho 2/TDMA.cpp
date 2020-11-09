#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include "alg_utilities.h"

double analytic_solution(double x, double t);
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& f);
void print_analytic();

//Dados:
const double alfa = 0.1;
const double L = 4.0;
const double T1 = 1.0;
const double T2 = 0.0;
const double t_max = 1.0;
const double beta = 0.5;
const double u = 0.5;
const double deltat = 0.05;
const double deltax = 0.2;
const double sigma = 0; //ou 1

int main(int argc, char **argv) {

	//Determinação dos coeficientes a serem utilizados:
	const auto C = u*deltat/deltax;
	const auto s = alfa*deltat/(deltax*deltax);
	std::cout << "C = " << C << ", s= " << s << std::endl;

	auto coef_a = 1+ 2*beta*(0.5*C*sigma + s);
	auto coef_b = beta*(0.5*C*(sigma-1) + s);
	auto coef_c = beta*(0.5*C*(sigma+1) + s);
	auto coef_d1 = (1-beta)*(0.5*C*(sigma+1) + s);
	auto coef_d2 = 1-2*(1-beta)*(0.5*C*sigma + s);
	auto coef_d3 = (1-beta)*(0.5*C*(sigma-1) + s);

	//Discretização em N intervalos:
	const size_t N = 11;

	//Inicialização dos vetores a serem utilizados:
	std::vector<double> a(N, coef_a);	//diagonal principal
	std::vector<double> b(N-1, coef_b);	//diagonal secundária superior
	std::vector<double> c(N-1, coef_c);	//diagonal secundária inferior

	std::vector<double> d(N, 0.0);
	std::vector<double> f(N, 0.0);	//vetor solução

	for (auto &e: a){
		std::cout << e << ' ';
	}
	std::cout << std::endl;
	for (auto &e: b){
		std::cout << e << ' ';
	}
	std::cout << std::endl;
	for (auto &e: c){
		std::cout << e << ' ';
	}
	std::cout << "\n\n";

	//Temperaturas iniciais fornecidas:
	f[0] = 1.0; f[10] = 0.0;

	std::cout << std::setprecision(6) << "T = ";
	for (size_t i=0; i<N; i++) {
		std::cout << f[i];
		if (i < N-1) {
			std:: cout << ' ';
		}
	}
	std::cout << std::endl;
	for (size_t i=1; i<N-1; i++) {
		d[i] = coef_d1*f[i-1] + coef_d2*f[i] + coef_d3*f[i+1] ;
	}

	tdma_solver(c, a, b, d, f);

	std::cout << "T = ";
	for (size_t i=0; i<N; i++) {
		std::cout << f[i];
		if (i < N-1) {
			std:: cout << ' ';
		}
	}
	std::cout << std::endl;
	print_analytic();
	return 0;
}

void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c,
std::vector<double>& d, std::vector<double>& f){
    auto N = d.size();
    const unsigned int n = N - 1;
    c[0] = c[0]/b[0];
    d[0] = d[0]/b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i > 0; i--){
        d[i] = d[i] - c[i]*d[i+1];
    }
}

void print_analytic(){
	std::fstream printer {"output_thomas_analytic.dat", std::ios::app};
	printer << "Solucao analitica.\n";
	printer << "Iteration x T\n";

	int i = 0;
	for (double x = -2;  x <= 2; x+= 0.2){
		double temp = analytic_solution(x, 1.0);
		printer << i << ' ' << x << ' ' << temp << "\n";
		i++;
	}
}

double analytic_solution(double x, double t){
	double sum = 0.0;
	for (int k = 1; k < 500; k++){
		sum += (1/(2*k - 1)) * std::exp((-alfa*(2*k - 1)*(2*k - 1)*NPI*NPI*t)/(L*L)) * std::sin((2*k - 1)*NPI*(x - u*t)/L);
	}
	return 0.5 - (2/NPI) * sum;
}
