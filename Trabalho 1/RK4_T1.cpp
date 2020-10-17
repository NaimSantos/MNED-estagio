#include <iostream>
#include <cmath>
#include <array>
#include <fstream>

//Dados:
const double Pr {10.0};
double X_0 {0.0}, Z_0 {0.5}, Y_0 {0.5};
const double b {8.0/3.0};
double Tao_0 {0.0};
const double Tao_f {50.0};
const double Tao_delt {0.05};
const double r{0.9};

double f1(double tao, double x, double y, double z);
double f2(double tao, double x, double y, double z);
double f3(double tao, double x, double y, double z);
void solveRK4(const double h, double t_i, const double t_f, std::array<double, 3>& S);

int main(){
	std::array<double, 3> S = {X_0, Y_0, Z_0};
	solveRK4(Tao_delt, Tao_0, Tao_f, S);
}
void solveRK4(const double h, double t_i, const double t_f, std::array<double, 3>& S){

	const auto n = std::round((t_f - t_i) / h);
	std::array<std::array<double, 5>, 3> k; //4 constantes do RK4 + 1 coluna para a soma dos ks

	std::fstream printer {"output_rk4.dat", std::ios::app};
	printer << "Solution via 4th Order Runge-Kutta Method.\n"; 
	printer << "t_inicial = " << t_i << ", total steps = " << n << ", t_final = " << t_f << '\n';
	printer << "Stept X Y Z\n";

	printer << "0 " << t_i << ' ' << S[0] << ' ' << S[1] << ' ' << S[2] << '\n';

	for (size_t i {1}; i <= n ; i++){
		
		k[0][0] = f1(t_i, S[0], S[1], S[2]);
		k[1][0] = f2(t_i, S[0], S[1], S[2]);
		k[2][0] = f3(t_i, S[0], S[1], S[2]);

		k[0][1] = f1(t_i + 0.5*h, S[0] + 0.5*k[0][0]*h, S[1] + 0.5*k[1][0]*h, S[2] + 0.5*k[2][0]*h);
		k[1][1] = f2(t_i + 0.5*h, S[0] + 0.5*k[0][0]*h, S[1] + 0.5*k[1][0]*h, S[2] + 0.5*k[2][0]*h);
		k[2][1] = f3(t_i + 0.5*h, S[0] + 0.5*k[0][0]*h, S[1] + 0.5*k[1][0]*h, S[2] + 0.5*k[2][0]*h);

		k[0][2] = f1(t_i + 0.5*h, S[0] + 0.5*k[0][1]*h, S[1] + 0.5*k[1][1]*h, S[2] + 0.5*k[2][1]*h);
		k[1][2] = f2(t_i + 0.5*h, S[0] + 0.5*k[0][1]*h, S[1] + 0.5*k[1][1]*h, S[2] + 0.5*k[2][1]*h);
		k[2][2] = f3(t_i + 0.5*h, S[0] + 0.5*k[0][1]*h, S[1] + 0.5*k[1][1]*h, S[2] + 0.5*k[2][1]*h);

		k[0][3] = f1(t_i + h, S[0] + k[0][2]*h, S[1] + k[1][2]*h, S[2] + k[2][2]*h);
		k[1][3] = f2(t_i + h, S[0] + k[0][2]*h, S[1] + k[1][2]*h, S[2] + k[2][2]*h);
		k[2][3] = f3(t_i + h, S[0] + k[0][2]*h, S[1] + k[1][2]*h, S[2] + k[2][2]*h);

		k[0][4] = (1.0/6.0)*(k[0][0] + 2*k[0][1] + 2*k[0][2] + k[0][3])*h; //Para a soma das constantes do RK4
		k[1][4] = (1.0/6.0)*(k[1][0] + 2*k[1][1] + 2*k[1][2] + k[1][3])*h; //Para a soma das constantes do RK4
		k[2][4] = (1.0/6.0)*(k[2][0] + 2*k[2][1] + 2*k[2][2] + k[2][3])*h; //Para a soma das constantes do RK4

		S[0] = S[0] + k[0][4];
		S[1] = S[1] + k[1][4];
		S[2] = S[2] + k[2][4];

		t_i = t_i + h;
		printer << i << ' ' << t_i << ' ' << S[0] << ' ' << S[1] << ' ' << S[2] << '\n';
	}
}
double f1(double tao, double x, double y, double z){
	return Pr*(y-x);
}
double f2(double tao, double x, double y, double z){
	return -x*z + r*x - y;
}
double f3(double tao, double x, double y, double z){
	return x*y - b*z;
}