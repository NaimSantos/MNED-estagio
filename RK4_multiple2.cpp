#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <iomanip>	//std::setprecision
//#include <functional>	//std::function

double f1(double x, double y, double z);
double f2(double x, double y, double z);
void solveRK4(const double h, double xi, const double xf, std::array<double, 2>& Ys);

int main(int argc, char* arg[]){
	double xi {0.0};
	const double xf {2.0};
	const double h {0.1};
	std::array<double, 2> Ys = {4.0, 6.0};
	
	solveRK4(h, xi, xf, Ys);

}
void solveRK4(const double h, double xi, const double xf, std::array<double, 2>& Ys){

	const auto n = std::round((xf - xi) / h);//numero de intervalos a estimar
	//double k[2][4+1] = {0.0}; //last collumn used for the sum of the constants
	std::array<std::array<double, 4+1>, 2> k;

	std::fstream printer {"output_rk4.csv", std::ios::app};
	printer << "Solution via 4th Order Runge-Kutta Method.\n"; 
	printer << "x_initial = " << xi << ", total steps = " << n << ", x_final = " << xf << '\n';
	printer << "step,t,y1,y2\n";
	
	for (size_t i {1}; i <= n ; i++){

		k[0][0] = f1(xi, Ys[0], Ys[1]);
		k[1][0] = f2(xi, Ys[0], Ys[1]);

		k[0][1] = f1(xi + 0.5*h, Ys[0] + 0.5*k[0][0]*h, Ys[1] + 0.5*k[1][0]*h);
		k[1][1] = f2(xi + 0.5*h, Ys[0] + 0.5*k[0][0]*h, Ys[1] + 0.5*k[1][0]*h);

		k[0][2] = f1(xi + 0.5*h, Ys[0] + 0.5*k[0][1]*h, Ys[1] + 0.5*k[1][1]*h);
		k[1][2] = f2(xi + 0.5*h, Ys[0] + 0.5*k[0][1]*h, Ys[1] + 0.5*k[1][1]*h);

		k[0][3] = f1(xi + h, Ys[0] + k[0][2]*h, Ys[1] + k[1][2]*h);
		k[1][3] = f2(xi + h, Ys[0] + k[0][2]*h, Ys[1] + k[1][2]*h);

		k[0][4] = (1.0/6.0)*(k[0][0] + 2*k[0][1] + 2*k[0][2] + k[0][3])*h;
		k[1][4] = (1.0/6.0)*(k[1][0] + 2*k[1][1] + 2*k[1][2] + k[1][3])*h;

		Ys[0] = Ys[0] + k[0][4];
		Ys[1] = Ys[1] + k[1][4];
		xi = xi + h;
		printer << i << ',' << xi << ',' << Ys[0] << ','  << Ys[1] << '\n';
		//std::cout << "x " << xi << ", y= " << Ys[0] << ", z= " << Ys[1] << std::endl;
	}
}
double f1(double x, double y, double z){
	return -0.5*y;
}
double f2(double x, double y, double z){
	return 4 - 0.3*z - 0.1*y;
}
