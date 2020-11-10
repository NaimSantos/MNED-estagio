#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>
#include <algorithm>	//std::fill
#include "alg_utilities.h"	//NPI

using std::vector;

double analytic_solution(double x, double t);
void tdma_solver(const vector<double>& a, const vector<double>& b, vector<double>& c, vector<double>& d, vector<double>& f, vector<double>& c_temp, vector<double>& d_temp);
void print_analytic();
void print_numerical(const vector<double>& f);

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
constexpr auto C = (u*deltat)/deltax;
constexpr auto s = (alfa*deltat)/(deltax*deltax);
constexpr auto coef_a = 1 + 2*beta*(0.5*C*sigma + s);
constexpr auto coef_b = beta*(0.5*C*(sigma-1) + s);
constexpr auto coef_c = beta*(0.5*C*(1+sigma) + s);
constexpr auto coef_d1 = (1-beta)*(0.5*C*(1+sigma) + s);
constexpr auto coef_d2 = 1 - 2*(1-beta)*(0.5*C*sigma + s);
constexpr auto coef_d3 = (1-beta)*(0.5*C*(sigma-1) + s);

int main(int argc, char **argv) {

	//Discretização em N intervalos:
	const size_t N = 11;
	//Inicialização dos vetores a serem utilizados:
	vector<double> a(N, coef_a);       //vetor da diagonal principal
	vector<double> b(N, coef_b);       //vetor da diagonal secundária superior
	vector<double> c(N, coef_c);       //vetor da diagonal secundária inferior
	vector<double> d(N, 0.0);          //vetor dos termos independentes
	vector<double> f(N, 0.0);          //vetor onde a solução será escrita
	vector<double> c_temp(N, 0.0);     //vetor auxiliar, para evitar que c seja sobreescrito
	vector<double> d_temp(N, 0.0);     //vetor auxiliar, para evitar que d seja sobreescrito

	b[N-1] = 0;	c[0] = 0; //Estes 2 vetores na verdade tem 1 elemento a menos que a diagonal principal
	d[0] = T1; d[N-1] = T2;	//Estes são os valores de temperatura conhecidos
	
	std::cout << std::setprecision(6) << "T = ";
	for (auto &e: f)
		std::cout << e << ' ';
	std::cout << std::endl;

	for (size_t i=1; i<N-1; i++) { //de d[1] ate d[10]
		d[i] = coef_d1*f[i-1] + coef_d2*f[i] + coef_d3*f[i+1] ;
	}

	tdma_solver(a, b, c, d, f, c_temp, d_temp);
	//tdma_solver(a1, b1, c1, d1, f1);

	std::cout << "T = ";
	for (auto &e: d)
		std::cout << e << ' ';

	//print_analytic();
	//print_numerical(f);
	return 0;
}
/*Parâmetros:
	a = diagonal inferior, b = diagonal principal, c = diagonal superior, d = termos independentes, f = vetor solução
*/
void tdma_solver(const vector<double>& a, const vector<double>& b, vector<double>& c, vector<double>& d, vector<double>& f, vector<double>& c_temp, vector<double>& d_temp){
	auto N = d.size() - 1;

	c[0] = c[0] / b[0];
	d[0] = d[0] / b[0];

	std::fill(c_temp.begin(), c_temp.end(), 0.0);	//Reseta os vetores temporarios, apenas por precaução
	std::fill(d_temp.begin(), d_temp.end(), 0.0);	//Reseta os vetores temporarios, apenas por precaução

	c_temp[0] = c[0] / b[0];
	d_temp[0] = d[0] / b[0];

	for (int i = 1; i < n; i++){
		c[i] = c[i] / (b[i] - a[i]*c[i-1]);
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
	}

	d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

	for (int i = n; i > 0; i--){
		d[i] = d[i] - c[i]*d[i+1];
	}
	//Escreve em f a solução:
	/*
	for (int i = n; i > 0; i--){
		f[i] = d[i] - c[i]*d[i+1];
	}
	*/
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
void print_numerical(const vector<double>& f){
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
