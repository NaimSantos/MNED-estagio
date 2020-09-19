#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <functional>
#include <iomanip>	//std::setprecision

double f1(double x, double y, double z);
double f2(double x, double y, double z);
void solveRK4(const double h, const double xi, const double xf, std::array<double, 2> Ys)

int main(int argc, char* arg[]){
	const double xi {0.0};
	const double xf {2.0};
	const double h {0.1};
	std::array<double, 2> Ys = {4.0, 6.0};
	
	solveRK4(h, xi, xf, Ys);
	std::cout << "y1:" << std::endl;
	std::cout << "z2:" << std::endl;
}
void solveRK4(const double h, const double xi, const double xf, std::array<double, 2> Ys){

	const auto n = std::round((xf - x1) / h);//numero de intervalos a estimar
	double ks[2][4+1] = {0.0}; //last collumn used for the sum of the slopes
	//std::array<std::array<double, 4+1>, 2> k;

	std::cout << "n= " << n << ", x0= " << x1 << ", y0= " << y1 << ", z0= " << y2 << ", xf=" << xf << std::endl;
	
	for (size_t i {1}; i <= n ; i++){

		k[0][0] = f1(x1, y1, y2);
		k[1][0] = f2(x1, y1, y2);

		k[0][1] = f1(x1 + 0.5*h, y1 + 0.5*k[0][0]*h, y2 + 0.5*k[1][0]*h);
		k[1][1] = f2(x1 + 0.5*h, y1 + 0.5*k[0][0]*h, y2 + 0.5*k[1][0]*h);

		k[0][2] = f1(x1 + 0.5*h, y1 + 0.5*k[0][1]*h, y2 + 0.5*k[1][1]*h);
		k[1][2] = f2(x1 + 0.5*h, y1 + 0.5*k[0][1]*h, y2 + 0.5*k[1][1]*h);

		k[0][3] = f1(x1 + h, y1 + k[0][2]*h, y2 + k[1][2]*h);
		k[1][3] = f2(x1 + h, y1 + k[0][2]*h, y2 + k[1][2]*h);

		k[0][4] = (1.0/6.0)*(k[0][0] + 2*k[0][1] + 2*k[0][2] + k[0][3])*h;
		l[1][4] = (1.0/6.0)*(k[1][0] + 2*k[1][1] + 2*k[1][2] + k[1][3])*h;

		y1 = y1 + k[0][4];
		y2 = y2 + l[1][4];
		x1 = x1 + h;
	}
}
double f1(double x, double y, double z){
	return -0.5*y;
}
double f2(double x, double y, double z){
	return 4 - 0.3*z - 0.1*y;
}