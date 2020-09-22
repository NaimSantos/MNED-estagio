#include <iostream>
#include <iomanip>	//std::setprecision
#include <cmath>

double solveRK4(double xi, double xf, double y0,  double step);
double solveRK3(double xi, double xf, double y0,  double step);
double solveRK5(double xi, double xf, double y0,  double step);
double solveEuler(double xi, double xf, double y0, double step);
double solveEuler_Mel(double xi, double xf, double y0, double step);
double solveEuler_Mod(double xi, double xf, double y0, double step);
double solveHuen_C(double xi, double xf, double y0, double step);
double solveRalston(double xi, double xf, double y0, double step);
double fxy(double x, double y);

int main(int argc, char* argv[]){
	
	//double x0{0.0}, x{0.2}, y0{1.0}, h{0.1};
	double x0{0.0}, x{2.0}, y0{0.5}, h{0.5};

	auto res1 = solveEuler(x0, x, y0, h);
	auto res2 = solveEuler_Mel(x0, x, y0, h);
	auto res3 = solveEuler_Mod(x0, x, y0, h);
	auto res4 = solveHuen_C(x0, x, y0, h);
	auto res5 = solveRalston(x0, x, y0, h);
	auto res6 = solveRK3(x0, x, y0, h);
	auto res7 = solveRK4(x0, x, y0, h);
	auto res8 = solveRK5(x0, x, y0, h);
	std::cout << "solveEuler: " << std::setprecision(12) << res1 << std::endl;
	std::cout << "solveEuler_Mel: " << res2 << std::endl;
	std::cout << "solveEuler_Mod: " << res3 << std::endl;
	std::cout << "solveHuen_C: " << res4 << std::endl;
	std::cout << "solveRalston: " << res5 << std::endl;
	std::cout << "solveRK3: " << res6 << std::endl;
	std::cout << "solveRK4: " << res7 << std::endl;
	std::cout << "solveRK5: " << res8 << std::endl;

}
double fxy(double x, double y){//dy/dx, a função conhecida
	//return 4*std::exp(0.8*x) - 0.5*y;
	//return - y + 1 - x;
	//return x*x*x*std::exp(-2*x) - 2*y;
	return y - x*x + 1; //exemplo de aula para rk4
}
double solveEuler(double xi, double xf, double y0, double h){
	auto n = std::round((xf - xi) / h) ;
	if (n <=0 )
		return NAN;

	for (int i = 1; i <= n; i++){
		y0 = y0 + fxy(xi, y0)*h;
		xi = xi + h;
	}
	return y0;
}
double solveEuler_Mel(double xi, double xf, double y0, double h){
	auto n = std::round((xf - xi) / h) ;
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;
	
	for (int i = 1; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + h, y + k1*h);

		y = y + (k1 + k2)*0.5*h;
		xi = xi + h;
	}
	return y;
}
double solveEuler_Mod(double xi,  double xf, double y0, double h){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.5*h, y + 0.5*k1*h);

		y = y + k2*h;
		xi = xi + h;
	}
	return y;
}
double solveHuen_C(double xi,  double xf, double y0, double h){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + h, y + k1*h);

		y = y + (0.5*k1 + 0.5*k2)*h;
		xi = xi + h;
	}
	return y;
}
double solveRalston(double xi,  double xf, double y0, double h){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0};
	double y = y0;

	for (size_t i {1}; i <= n; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.75*h, y + 0.75*k1*h);

		y = y + ((1.0/3.0)*k1 + (2.0/3.0)*k2)*h;
		xi = xi + h;
	}
	return y;
}
double solveRK3(double xi,  double xf, double y0, double h){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0}, k3{0};
	double y = y0;

	for (size_t i {1}; i <=n ; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.5*h, y + 0.5*k1*h);
		k3 = fxy(xi + h, y - k1*h + 2*k2*h);

		y = y + (1.0/6.0)*(k1 + 4*k2 + k3)*h;
		xi = xi + h;
	}
	return y;
}
double solveRK4(double xi,  double xf, double y0, double h){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0}, k3{0}, k4{0};
	double y = y0;

	for (size_t i {1}; i <= n ; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.5*h, y + 0.5*k1*h);
		k3 = fxy(xi + 0.5*h, y + 0.5*k2*h);
		k4 = fxy(xi + h, y + k3*h);

		y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)*h;
		xi = xi + h;
	}
	return y;
}
double solveRK5(double xi,  double xf, double y0, double h){
	auto n = std::round((xf - xi) / h);
	if (n <=0 )
		return NAN;

	double k1{0}, k2{0}, k3{0}, k4{0}, k5{0}, k6{0};
	double y = y0;

	for (size_t i {1}; i <=n ; i++){
		k1 = fxy(xi, y);
		k2 = fxy(xi + 0.25*h, y + 0.25*k1*h);
		k3 = fxy(xi + 0.25*h, y + 0.125*k1*h + 0.125*k2*h);
		k4 = fxy(xi + 0.5*h, y - 0.5*k2*h + k3*h);
		k5 = fxy(xi + 0.75*h, y + (3.0/16.0)*k1*h + (9.0/16.0)*k4*h);
		k6 = fxy(xi + h, y - (3.0/7.0)*k1*h + (2.0/7.0)*k2*h + (12.0/7.0)*k3*h - (12.0/7.0)*k4*h + (8.0/7.0)*k5*h);

		y = y + (1.0/90.0)*(7*k1 + 32*k2 + 12*k4 + 32*k5 + 7*k6)*h;
		xi = xi + h;
	}
	return y;
}