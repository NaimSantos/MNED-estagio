#include <iostream>
#include <cmath>	//para funções trigonometricas, logaritmicas e exponencial
#include <array>
#include <fstream> //std::fstream
#include <iomanip>	//std::setprecision
//#include <functional>	//std::function

double f1(double x, double y, double z);
double f2(double x, double y, double z);
template <size_t N>	void solveRK4(const double h, double xi, const double xf, std::array<double, N>& Y);

int main(int argc, char* arg[]){

	double xi {0.0};
	const double xf {0.3};
	const double h {0.1};
	
	std::array<double, 2> Y = {0.0, 0.0};
	
	solveRK4<Y.size()>(h, xi, xf, Y); //Y.size() --> size_t
}
template <size_t N>//qualquer N 
void solveRK4(const double h, double xi, const double xf, std::array<double, N>& Y){

	const auto n = std::round((xf - xi) / h);//numero de intervalos a estimar
	//double k[2][4+1] = {0.0}; //ultima coluna para a soma dos ks
	std::array<std::array<double, 4+1>, N> k;

	std::fstream printer {"output_rk4.csv", std::ios::app};
	printer << "Solution via 4th Order Runge-Kutta Method.\n"; 
	printer << "x_initial = " << xi << ", total steps = " << n << ", x_final = " << xf << '\n';
	printer << "step,t,y1,y2\n";
	
	for (size_t i {1}; i <= n ; i++){

		k[0][0] = f1(xi, Y[0], Y[1]);
		k[1][0] = f2(xi, Y[0], Y[1]);

		k[0][1] = f1(xi + 0.5*h, Y[0] + 0.5*k[0][0]*h, Y[1] + 0.5*k[1][0]*h);
		k[1][1] = f2(xi + 0.5*h, Y[0] + 0.5*k[0][0]*h, Y[1] + 0.5*k[1][0]*h);

		k[0][2] = f1(xi + 0.5*h, Y[0] + 0.5*k[0][1]*h, Y[1] + 0.5*k[1][1]*h);
		k[1][2] = f2(xi + 0.5*h, Y[0] + 0.5*k[0][1]*h, Y[1] + 0.5*k[1][1]*h);

		k[0][3] = f1(xi + h, Y[0] + k[0][2]*h, Y[1] + k[1][2]*h);
		k[1][3] = f2(xi + h, Y[0] + k[0][2]*h, Y[1] + k[1][2]*h);

		k[0][4] = (1.0/6.0)*(k[0][0] + 2*k[0][1] + 2*k[0][2] + k[0][3])*h; //Para a soma das constantes
		k[1][4] = (1.0/6.0)*(k[1][0] + 2*k[1][1] + 2*k[1][2] + k[1][3])*h; //Para a soma das constantes

		Y[0] = Y[0] + k[0][4];
		Y[1] = Y[1] + k[1][4];
		xi = xi + h;
		printer << i << ',' << xi << ',' << Y[0] << ','  << Y[1] << '\n';
		std::cout << "x= " << std::setprecision(10) << xi << ", y= " << Y[0] << ", z= " << Y[1] << std::endl;
	}
}
double f1(double x, double y, double z){
	//return -0.5*y;
	return z;
}
double f2(double x, double y, double z){
	//return 4 - 0.3*z - 0.1*y;
	return 1 - z - y;
}