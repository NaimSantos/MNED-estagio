#include <iostream>
#include <cmath>

double solveBisect(double a, double b);
double ftarget(double x);

int main (int argc, char* argv[]){
	double a0{0}, b0{1.0};
	std::cout << "Informe valores para o intervalo de busca: ";
	std::cin >> a0 >> b0;
	auto res = solveBisect(a0, b0);
	if (std::isnan(res))
		std::cout << "Nao foi encontrada uma raiz" << std::endl;
}
double ftarget(double x){
	return (x * x * x) - (3 * x * x) + 2;
}
double solveBisect(double a0, double b0){
	if (ftarget(a0) * ftarget(b0) >= 0){
		return NAN;
	}
	double c0 {a0};
	double eps {0.00001};
	size_t i{1};
	while ((b0-a0) >= eps || i<100){
		c0 = (a0 + b0) / 2;
		if (ftarget(c0) == 0)
			break;
		else
			(ftarget(c0)*ftarget(a0) < 0) ? (b0 = c0) : (a0 = c0);
		i++;
	}
	if(i>=100){
		return NAN;
	}
	std::cout << "\nRaiz " << c0 <<" encontrada na interacao " << i-1 << std::endl;
	return c0;
}