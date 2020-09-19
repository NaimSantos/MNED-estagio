#include <iostream>
#include <cmath>
#include <array>
#include <iomanip>	//std::setprecision
#include <utility> //std::pair

double fxy(double x, double y, double z);
double f1(double x, double y, double z);
double f2(double x, double y, double z);
void solveRK4(double h, double xf, std::pair<double, double>& xy1, std::pair<double, double>& xy2);

int main(int argc, char* arg[]){
	std::pair<double, double> xy1 = {0.0, 4.0};
	std::pair<double, double> xy2 = {xy1.first, 6.0};

	double xf {2.0};
	double h {0.1};
	
	solveRK4(h, xf, xy1, xy2);
	std::cout << "x:" << xy1.first << std::endl;
	std::cout << "y:" << xy1.second << std::endl;
	std::cout << "z:" << xy2.second << std::endl;
}
void solveRK4(double h, double xf, std::pair<double, double>& xy1, std::pair<double, double>& xy2){
	auto& x1 = xy1.first;
	auto& y1 = xy1.second;
	//auto& x2 = xy2.first;
	auto& y2 = xy2.second;
	auto n = std::round((xf - x1) / h);//numero de intervalos a estimar

	double k1{0}, k2{0}, k3{0}, k4{0};
	double l1{0}, l2{0}, l3{0}, l4{0};
	double k {0.0}, l {0.0};
	std::cout << "n= " << n << ", x0= " << x1 << ", y0= " << y1 << ", z0= " << y2 << ", xf=" << xf << std::endl;
	
	for (size_t i {1}; i <= n ; i++){

		k1 = f1(x1, y1, y2);
		l1 = f2(x1, y1, y2);

		k2 = f1(x1 + 0.5*h, y1 + 0.5*k1*h, y2 + 0.5*l1*h);
		l2 = f2(x1 + 0.5*h, y1 + 0.5*k1*h, y2 + 0.5*l1*h);

		k3 = f1(x1 + 0.5*h, y1 + 0.5*k2*h, y2 + 0.5*l2*h);
		l3 = f2(x1 + 0.5*h, y1 + 0.5*k2*h, y2 + 0.5*l2*h);

		k4 = f1(x1 + h, y1 + k3*h, y2 + l3*h);
		l4 = f2(x1 + h, y1 + k3*h, y2 + l3*h);

		k = (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)*h;
		l = (1.0/6.0)*(l1 + 2*l2 + 2*l3 + l4)*h;

		y1 = y1 + k;
		y2 = y2 + l;
		x1 = x1 + h;
	}
}
double f1(double x, double y, double z){
	return -0.5*y;
}
double f2(double x, double y, double z){
	return 4 - 0.3*z - 0.1*y;
}