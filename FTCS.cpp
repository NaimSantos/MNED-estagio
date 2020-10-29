#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>

#include "time_utilities.h"

int i4_modp(int i, int j);
int i4_wrap(int ival, int ilo, int ihi);
double* initial_condition(int nx, double x[]);
double* r8vec_linspace_new(int n, double a, double b);
void writegnuplot(const std::string& filename);

int main(){
	print_current_time();
	auto start = capturetime();


	std::string data_filename = "dados_adveccao.txt";
	std::ofstream data_output;

	double a, b, c;
	double dt, dx;
	int jm1, jp1, nx, nt, nt_step, t;
	double *u;
	double *unew;
	double *x;

	//Definição de variáveis:
	nx = 101;
	dx = static_cast<double>(1.0 / (nx-1));
	a = 0.0;
	b = 1.0;
	x = r8vec_linspace_new(nx, a, b);
	nt = 1000;
	dt = static_cast<double>(1.0 / nt);
	c = 1.0;

	//Define as condições iniciais:
	u = initial_condition(nx, x);

	//Abre o arquivo e escreve as soluções enquanto elas são computadas:.
	data_output.open(data_filename);

	t = 0.0;
	data_output << "  " << x[0] << "  " << t << "  " << u[0] << "\n";
	for (int j = 0; j < nx; j++ ){
		data_output << "  " << x[j] << "  " << t << "  " << u[j] << "\n";
	}
	data_output << "\n";

	nt_step = 100;
	unew = new double[nx];

	for (int i = 0; i < nt; i++ ){
		for (int j = 0; j < nx; j++ ){
			jm1 = i4_wrap (j-1, 0, nx-1);
			jp1 = i4_wrap (j+1, 0, nx-1);
			unew[j] = u[j] - c * dt / dx / 2.0 * (u[jp1] - u[jm1]);
		}
		for (int j = 0; j < nx; j++){
			u[j] = unew[j];
		}
		if (i == nt_step-1){
			t =  i*dt;
			for (int j = 0; j < nx; j++){
				data_output << "  " << x[j] << "  " << t << "  " << u[j] << "\n";
			}
			data_output << "\n";
			nt_step = nt_step + 100;
		}
	}

	//Fechar o arquivo quando acabar de computar:
	data_output.close();
	std::cout << "\n\tResultados exportados para \"" << data_filename << "\"\n";

	//Gerar um arquivo com comandos gnuplot:
	std::string gnuplotfile = "comandos_adveccao.plt";
	writegnuplot(gnuplotfile);
	std::cout << "\tComandos do Gnuplot exportados para \"" << gnuplotfile << "\"\n\n";

	//Desalocar as matrizes:
	delete [] u;
	delete [] unew;
	delete [] x;

	//Encerrar:
	auto end = capturetime();
	auto total = get_elapsed_time(start, end);
	std::cout << "\nTotal execution time: " << total * 1e-06 << " s"<<std::endl;
	return 0;
}

int i4_modp (int i, int j){

	if (j == 0){
		std::cerr << "\n";
		std::cerr << "I4_MODP - Erro!\n";
		std::cerr << "I4_MODP (I, J) chamado with J = " << j << "\n";
		std::exit(1);
	}

	int value = i % j;

	if (value < 0){
		value = value + std::abs(j);
	}

	return value;
}

int i4_wrap (int ival, int ilo, int ihi){
	int jhi, jlo, value, wide;

	if ( ilo<=ihi ){
		jlo = ilo;
		jhi = ihi;
	}
	else {
		jlo = ihi;
		jhi = ilo;
	}

	wide = jhi + 1 - jlo;

	(wide==1) ? (value = jlo) : (value = jlo + i4_modp (ival-jlo, wide));

	return value;
}

double* initial_condition(int nx, double x[]){

	double *u = new double[nx];

	for (int i = 0; i < nx; i++ ){
		if ( 0.4<=x[i] && x[i]<=0.6 ){
			u[i] = std::pow( 10.0*x[i] - 4.0, 2) * std::pow(6.0 - 10.0*x[i], 2);
		}
		else {
			u[i] = 0.0;
		}
	}

	return u;
}

double* r8vec_linspace_new(int n, double a_first, double a_last){
	double *a = new double[n];

	if (n==1){
		a[0] = ( a_first + a_last ) / 2.0;
	}
	else {
		for (int i = 0; i < n; i++ ){
			a[i] = (( n-1-i ) * a_first + (i) * a_last) / ( n-1);
		}
	}
	return a;
}

void writegnuplot(const std::string& filename){
	std::ofstream fwritter;

	fwritter.open(filename);

	fwritter << "set term png\n";
	fwritter << "set output 'AdvectionGraph.png'\n";
	fwritter << "set grid\n";
	fwritter << "set style data lines\n";
	fwritter << "unset key\n";
	fwritter << "set xlabel 'X'\n";
	fwritter << "set ylabel 'Time'\n";
	fwritter << "splot '" << filename << "' using 1:2:3 with lines\n";
	fwritter << "quit\n";

	fwritter.close();
}