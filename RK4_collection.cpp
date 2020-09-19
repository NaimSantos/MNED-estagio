#include <iostream>
#include <iomanip>	//std::setprecision
#include <cmath>

double solveRK4(double xi, double xf, double y0,  double step);
double solveRK3(double xi, double xf, double y0,  double step);
double solveRK5(double xi, double xf, double y0,  double step);
double solveEuler(double xi, double xf, double y0, double step);
double solveEuler_M(double xi, double xf, double y0, double step);
double solvePonto_Med(double xi, double xf, double y0, double step);
double solveHuen_C(double xi, double xf, double y0, double step);
double solveRalston(double xi, double xf, double y0, double step);
double fxy(double x, double y);

int main(int argc, char* argv[]){

	//double x0{0}, x{0}, y0{0}, h{0};
	//std::cout << "Informe valores iniciais (x0), final (x), o passo (h) e o ponto de busca (y0):\n";
	//std::cin>> x0 >> x >> h >> y0 ;
	double x0{0}, x{4}, y0{2}, h{0.01};
	auto res1 = solveEuler(x0, x, y0, h);
	auto res3 = solveRK3(x0, x, y0, h);
	auto res4 = solveRK4(x0, x, y0, h);
	auto res5 = solveRK5(x0, x, y0, h);
	std::cout << "Euler: " << std::setprecision(10) << res1 << std::endl;
	std::cout << "solveRK3: " << res3 << std::endl;
	std::cout << "solveRK4: " << res4 << std::endl;
	std::cout << "solveRK5: " << res5 << std::endl;
}
double fxy(double x, double y){//dy/dx, a função conhecida
	return 4*std::exp(0.8*x) - 0.5*y;
}
double solveEuler(double xi, double xf, double y0, double h){
	while (xi < xf) {
		y0 = y0 + fxy(xi, y0)*h;
		//std::cout << "x= " << xi << " , y= " << y0 <<std::endl;
		xi = xi + h;
	}
	return y0;
}
double solveEuler_M(double xi, double xf, double y0, double h){
	auto n = ((xf - xi) / h) ;
	for (int i = 0; i <= n; i++){
		y0 = y0 + fxy(xi + i*h, y0 + i*h)*h;
	}
	return y0;
}
double solvePonto_Med(double x0,  double x, double y0, double h){
	auto n = static_cast<unsigned int>((x - x0) / h);
	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(x0, y);
		k2 = fxy(x0 + 0.5*h, y + 0.5*k1*h);

		y = y + k2*h;
		x0 = x0 + h;
	}
	return y;
}
double solveHuen_C(double x0,  double x, double y0, double h){
	auto n = static_cast<unsigned int>((x - x0) / h);
	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(x0, y);
		k2 = fxy(x0 + h, y + k1*h);

		y = y + (0.5*k1 + 0.5*k2)*h;
		x0 = x0 + h;
	}
	return y;
}
double solveRalston(double x0,  double x, double y0, double h){
	auto n = static_cast<unsigned int>((x - x0) / h);
	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(x0, y);
		k2 = fxy(x0 + 0.75*h, y + 0.75*k1*h);

		y = y + ((1.0/3.0)*k1 + (2.0/3.0)*k2)*h;
		x0 = x0 + h;
	}
	return y;
}
double solveRK3(double x0,  double x, double y0, double h){
	auto n = static_cast<unsigned int>((x - x0) / h);
	double k1{0}, k2{0}, k3{0};
	double y = y0;

	for (size_t i {1}; i <=n ; i++){
		k1 = fxy(x0, y);
		k2 = fxy(x0 + 0.5*h, y + 0.5*k1*h);
		k3 = fxy(x0 + h, y - k1*h + 2*k2*h);

		y = y + (1.0/6.0)*(k1 + 4*k2 + k3)*h;
		x0 = x0 + h;
	}
	return y;
}
double solveRK4(double x0,  double x, double y0, double h){
	auto n = static_cast<unsigned int>((x - x0) / h);
	double k1{0}, k2{0}, k3{0}, k4{0};
	double y = y0;

	for (size_t i {1}; i <= n ; i++){
		k1 = fxy(x0, y);
		k2 = fxy(x0 + 0.5*h, y + 0.5*k1*h);
		k3 = fxy(x0 + 0.5*h, y + 0.5*k2*h);
		k4 = fxy(x0 + h, y + k3*h);

		y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)*h;
		x0 = x0 + h;
	}
	return y;
}
double solveRK5(double x0,  double x, double y0, double h){
	auto n = static_cast<unsigned int>((x - x0) / h);
	double k1{0}, k2{0}, k3{0}, k4{0}, k5{0}, k6{0};
	double y = y0;

	for (size_t i {1}; i <=n ; i++){
		k1 = fxy(x0, y);
		k2 = fxy(x0 + 0.25*h, y + 0.25*k1*h);
		k3 = fxy(x0 + 0.25*h, y + 0.125*k1*h + 0.125*k2*h);
		k4 = fxy(x0 + 0.5*h, y - 0.5*k2*h + k3*h);
		k5 = fxy(x0 + 0.75*h, y + (3.0/16.0)*k1*h + (9.0/16.0)*k4*h);
		k6 = fxy(x0 + h, y - (3.0/7.0)*k1*h + (2.0/7.0)*k2*h + (12.0/7.0)*k3*h - (12.0/7.0)*k4*h + (8.0/7.0)*k5*h);

		y = y + (1.0/90.0)*(7*k1 + 32*k2 + 12*k4 + 32*k5 + 7*k6)*h;
		x0 = x0 + h;
	}
	return y;
}