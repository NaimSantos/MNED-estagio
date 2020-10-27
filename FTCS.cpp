#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#include <cstdlib>
#include <ctime>
#include <cstring>


int i4_modp(int i, int j);
int i4_wrap(int ival, int ilo, int ihi);
double *initial_condition(int nx, double x[]);
double *r8vec_linspace_new(int n, double a, double b);
void time_capture();

int main (){
	double a, b, c;
	std::string command_filename = "comandos_adveccao.txt";
	std::ofstream command_unit;
	std::string data_filename = "dados_adveccao.txt";
	std::ofstream data_unit;
	double dt, dx;
	int jm1, jp1, nx, nt, nt_step, t;
	double *u;
	double *unew;
	double *x;

	time_capture();

	//Definição de variáveis:
	nx = 101;
	dx = static_cast<double>(1.0 / (nx-1));
	a = 0.0, 
	b = 1.0;
	x = r8vec_linspace_new(nx, a, b);
	nt = 1000;
	dt = static_cast<double>(1.0 / nt);
	c = 1.0;

	//Define as condições iniciais:
	u = initial_condition(nx, x);

	//Abre o arquivo e escreve as soluções enquanto elas são computadas:.
	data_unit.open(data_filename.c_str ());

	t = 0.0;
	data_unit << "  " << x[0] << "  " << t << "  " << u[0] << "\n";
	for (int j = 0; j < nx; j++ ){
		data_unit << "  " << x[j] << "  " << t << "  " << u[j] << "\n";
	}
	data_unit << "\n";

	nt_step = 100;

	std::cout << "\n";
	std::cout << "  Numero de nos NX = " << nx << "\n";
	std::cout << "  Numero de intervalos de tempo NT = " << nt << "\n";
	std::cout << "  Velocidade constant c = " << c << "\n";

	unew = new double[nx];

	for (int i = 0; i < nt; i++ ){
		for (int j = 0; j < nx; j++ ){
			jm1 = i4_wrap (j-1, 0, nx-1);
			jp1 = i4_wrap (j+1, 0, nx-1 );
			unew[j] = u[j] - c * dt / dx / 2.0 * (u[jp1] - u[jm1]);
		}
		for (int j = 0; j < nx; j++){
			u[j] = unew[j];
		}
		if (i == nt_step-1){
			t =  i*dt;
			for (int j = 0; j < nx; j++){
				data_unit << "  " << x[j] << "  " << t << "  " << u[j] << "\n";
			}
			data_unit << "\n";
			nt_step = nt_step + 100;
		}
	}

	//Fechar o arquivo quando acabar de computar
	data_unit.close();

	std::cout << "\n  Plotar valores em um arquivo \"" << data_filename << "\"\n";

	//Escreve a arquivo do gnuplot:

	command_unit.open(command_filename.c_str());

	command_unit << "set term png\n";
	command_unit << "set output 'advection.png'\n";
	command_unit << "set grid\n";
	command_unit << "set style data lines\n";
	command_unit << "unset key\n";
	command_unit << "set xlabel '<---X--->'\n";
	command_unit << "set ylabel '<---Time--->'\n";
	command_unit << "splot '" << data_filename << "' using 1:2:3 with lines\n";
	command_unit << "quit\n";

	command_unit.close ();

	std::cout << "  Comandos do Gnuplot foram escritos em um arquivo \"" << command_filename << "\"\n";

	//Desalocar as matrizes:

	delete [] u;
	delete [] unew;
	delete [] x;

	std::cout << "Fim de execucao normal.\n\n";
	time_capture();

	return 0;
}

int i4_modp (int i, int j){
	int value;

	if (j == 0){
		std::cerr << "\n";
		std::cerr << "I4_MODP - Erro!\n";
		std::cerr << "I4_MODP (I, J) chamado with J = " << j << "\n";
		exit (1);
	}

	value = i % j;

	if (value < 0){
		value = value + abs ( j );
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

double *initial_condition(int nx, double x[]){

	double *u = new double[nx];

	for (int i = 0; i < nx; i++ ){
		if ( 0.4 <= x[i] && x[i] <= 0.6 ){
			u[i] = std::pow( 10.0 * x[i] - 4.0, 2 ) * std::pow ( 6.0 - 10.0 * x[i], 2 );
		}
		else {
			u[i] = 0.0;
		}
	}

	return u;
}

double *r8vec_linspace_new(int n, double a_first, double a_last){
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

void time_capture(){
	const int time_size {40};
	static char time_buffer[time_size];
	const struct std::tm *tm_ptr;
	std::time_t now;

	now = std::time(NULL);
	tm_ptr = std::localtime (&now);

	std::strftime (time_buffer, time_size, "%d %B %Y %I:%M:%S %p", tm_ptr);

	std::cout << time_buffer << "\n";
}
