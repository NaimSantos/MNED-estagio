#include <iostream>
#include <cmath>

double solveNewton(double x);
double ftarget(double x);
double ftargetdiff(double x);

int main(int argc, char* argv[]){

	double x0 {0};
	std::cout << "Informe uma estimativa inicial para a raiz: ";
	std::cin >> x0 ;

	auto res = solveNewton(x0);
	std::cout << "A estimativa encontrada foi: " << res << std::endl;
}
double ftarget(double x){
	return std::pow(x, 5) + (3 * std::pow(x, 2)) - 1;
}
double ftargetdiff(double x){
	return 5 * std::pow(x, 4) + (6 * x);
}
double solveNewton(double x){
	auto h = ftarget(x) / ftargetdiff(x);
	unsigned int i {0};
	while (std::abs(h) >= 0.00001 && i<100){
		h = ftarget(x) / ftargetdiff(x);
		x = x- h;
		i++;
	}
	return x;
}