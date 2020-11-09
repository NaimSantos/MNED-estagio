#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>	//std::fill
#include "alg_utilities.h"	//NPI

double analytic_solution(double x, double t);
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& f, std::vector<double>& c_temp, std::vector<double>& d_temp);
void print_analytic();
void print_numerical(const std::vector<double>& f);

//Dados:
constexpr double alfa = 0.1;
constexpr double L = 4.0;
constexpr double T1 = 1.0;
constexpr double T2 = 0.0;
constexpr double t_max = 1.0;
constexpr double beta = 0.5;
constexpr double u = 0.5;
constexpr double deltat = 0.05;
constexpr double deltax = 0.2;
constexpr double sigma = 0; //0 ou 1

//Determinação dos coeficientes a serem utilizados:
constexpr auto C = u*deltat/deltax;
constexpr auto s = alfa*deltat/(deltax*deltax);
constexpr auto coef_a = 1+ 2*beta*(0.5*C*sigma + s);
constexpr auto coef_b = beta*(0.5*C*(sigma-1) + s);
constexpr auto coef_c = beta*(0.5*C*(sigma+1) + s);
constexpr auto coef_d1 = (1-beta)*(0.5*C*(sigma+1) + s);
constexpr auto coef_d2 = 1-2*(1-beta)*(0.5*C*sigma + s);
constexpr auto coef_d3 = (1-beta)*(0.5*C*(sigma-1) + s);

int main(int argc, char **argv) {

	//Discretização em N intervalos:
	const size_t N = 11;
	//Inicialização dos vetores a serem utilizados:
	std::vector<double> a(N, coef_a);	//diagonal principal
	std::vector<double> b(N-1, coef_b);	//diagonal secundária superior
	std::vector<double> c(N-1, coef_c);	//diagonal secundária inferior

	std::vector<double> d(N, 0.0);
	std::vector<double> f(N, 0.0);	//vetor solução
	std::vector<double> c_temp(N, 0.0);	//temporarios, para evitar que os vetores sejam sobreescritos
	std::vector<double> d_temp(N, 0.0);	//temporarios, para evitar que os vetores sejam sobreescritos

	//Temperaturas iniciais fornecidas:
	f[0] = 1.0;	d[0] = 1.0;

	std::cout << std::setprecision(6) << "T = ";
	for (auto &e: f){
		std::cout << e << ' ';
	}
	std::cout << std::endl;

	for (size_t i=1; i<N-1; i++) { //de d[1] ate d[10]
		d[i] = coef_d1*f[i-1] + coef_d2*f[i] + coef_d3*f[i+1] ;
	}

	tdma_solver(c, a, b, d, f, c_temp, d_temp);

	std::cout << "T = ";
	for (auto &e: f){
		std::cout << e << ' ';
	}

	print_analytic();
	print_numerical(f);
	return 0;
}
/*Parâmetros:
	a = diagonal inferior,
	b = diagonal principal
	c = diagonal superior
	d = termos independentes
	f = vetor solução
	c_temp e d_temp = vetores temporários para evitar que os vetores principais sejam sobreescrios.
*/
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& f, std::vector<double>& c_temp, std::vector<double>& d_temp){
	auto N = d.size();

	std::fill(c_temp.begin(), c_temp.end(), 0.0);
	std::fill(d_temp.begin(), d_temp.end(), 0.0);

	c_temp[0] = c[0] / b[0];
	d_temp[0] = d[0] / b[0];

	for (int i = 1; i < N; i++){
		double m = 1.0 / (b[i] - a[i] * c_temp[i-1]);
		c_temp[i] = c[i] * m;
		d_temp[i] = (d[i] - a[i] * d_temp[i-1]) * m;
	}
	for (int i = N-1; i > 0; i--){
		f[i] = d_temp[i] - c_temp[i] * d[i+1];
	}
}

void print_analytic(){
	std::fstream printer {"output_thomas_analytic.dat", std::ios::app};
	printer << "Solucao analitica.\n";
	printer << "Iteration x T\n";

	int i = 0;
	for (double x = -2;  x <= 2; x+= 0.4){
		double temp = analytic_solution(x, 1.0);
		printer << i << ' ' << x << ' ' << temp << "\n";
		i++;
	}
}
void print_numerical(const std::vector<double>& f){
	std::fstream printer {"output_thomas_crank_nicolson.dat", std::ios::app};
	printer << "Solucao Via Crank-Nicolson.\n";
	printer << "Iteration x T\n";

	int i = 0;
	for (double x = -2;  x <= 2; x+= 0.4){
		double temp = analytic_solution(x, 1.0);
		printer << i << ' ' << x << ' ' << f[i] << "\n";
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


/*
	//Exibir antes da execução:
	std::cout << "A: ";
	for (auto &e: a){
		std::cout << e << ' ';
	}
	std::cout << "\nB: ";
	for (auto &e: b){
		std::cout << e << ' ';
	}
	std::cout << "\nC: ";
	for (auto &e: c){
		std::cout << e << ' ';
	}
	std::cout << "\nD: ";
	for (auto &e: d){
		std::cout << e << ' ';
	}
	std::cout << "\n\n";

	//Exibir após a execução:
	std::cout << "\nA: ";
	for (auto &e: a){
		std::cout << e << ' ';
	}
	std::cout << "\nB: ";
	for (auto &e: b){
		std::cout << e << ' ';
	}
	std::cout << "\nC: ";
	for (auto &e: c){
		std::cout << e << ' ';
	}
	std::cout << "\nD: ";
	for (auto &e: d){
		std::cout << e << ' ';
	}
	std::cout << "\n\n";

*/