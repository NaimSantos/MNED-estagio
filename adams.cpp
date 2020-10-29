#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <fstream>
#include <iomanip>	//std::setprecision

double f(double x, double y);
void solveRK4(const double h, double xi, double xf, std::vector<double>& Y);
void solveAdams(const double h, double xi, const double xf, std::vector<double>& Y);

int main(int argc, char* arg[]){
	double xi {0.0};
	const double xf {0.8};
	const double h {0.2};

	std::vector<double> Y = {0.0};
	solveAdams(h, xi, xf, Y);
}
void solveRK4(const double h, double xi, double xf, std::vector<double>& Y){
	auto n = std::round((xf - xi) / h);
	double k1{0}, k2{0}, k3{0}, k4{0};
	double y = Y[0];

	for (size_t i {1}; i <= n ; i++){
		k1 = f(xi, y);
		k2 = f(xi + 0.5*h, y + 0.5*k1*h);
		k3 = f(xi + 0.5*h, y + 0.5*k2*h);
		k4 = f(xi + h, y + k3*h);

		y = y + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)*h;
		xi = xi + h;
		Y.push_back(y);
	}
}
void solveAdams(const double h, double xi, const double xf, std::vector<double>& Y){

	solveRK4(h, xi, xf-h, Y); //Runge kutta determina os valores não fornecidos
	for (auto& e : Y)
		std::cout << e << ' ';
	std::cout << std::endl;

	const auto n = std::round((xf - xi) / h);
	std::cout << "Numero de iteracoes estimadas: " << n << std::endl;

	double milne_p {0.0}, milne_c {0.0};
	auto x = xi;
	for (size_t i = 4; i <= n; ++i){

		auto f4 = f(x+i*h, Y[i]);
		auto f3 = f(x+(i-1)*h, Y[i-1]);
		auto f2 = f(x+(i-2)*h, Y[i-2]);
		auto f1 = f(x+(i-3)*h, Y[i-3]);

		milne_p = Y[i] + ((55.0/24.0)*f4 - (59.0/24.0)*f3 + (37.0/24.0)*f2 - (9.0/24.0)*f1)*h;
		milne_c = Y[i] + ((9.0/24.0)*f(x+i*h, milne_p) + (19.0/24.0)*f(x+(i)*h, Y[i]) - (5.0/24.0)*f(x+(i-1)*h, Y[i-1]) + (1.0/24.0)*f(x+(i-2)*h, Y[i-2]))*h;

		Y.push_back(milne_c);
		std::cout << "i:" << i << ", yp: " << milne_p << ", yc: " << milne_c << std::endl;
	}
}

//Funções avaliadas
double f(double x, double y){
	return 1 + y*y;
}
