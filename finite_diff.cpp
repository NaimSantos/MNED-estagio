#include <iostream>
#include <cmath>
#include <fstream>

double f(double x);

double fdiff_a_1st(double x, double h);
double fdiff_a_1st_v2(double x, double h);
double fdiff_a_2nd(double x, double h);
double fdiff_a_2nd_v2(double x, double h);

double fdiff_r_1st(double x, double h);
double fdiff_r_1st_v2(double x, double h);
double fdiff_r_2nd(double x, double h);
double fdiff_r_2nd_v2(double x, double h);

double fdiff_c_1st(double x, double h);
double fdiff_c_1st_v2(double x, double h);
double fdiff_c_2nd(double x, double h);
double fdiff_c_2nd_v2(double x, double h);

int main(int argc, char* arg[]){
	double h = 0.1;
	double x = 2.0;

	auto adiff1 = fdiff_a_1st(x, h);
	auto adiff2 = fdiff_a_1st_v2(x, h);
	auto rdiff1 = fdiff_r_1st(x, h);
	auto rdiff2 = fdiff_r_1st_v2(x, h);
	auto cdiff1 = fdiff_c_1st(x, h);
	auto cdiff2 = fdiff_c_1st_v2(x, h);

	std::cout <<  "Resultado por diferencas finitas avancada 1 : " << adiff1 << std::endl;
	std::cout <<  "Resultado por diferencas finitas avancada 2 : " << adiff2 << std::endl;
	std::cout <<  "Resultado por diferencas finitas recuada 1 : " << rdiff1 << std::endl;
	std::cout <<  "Resultado por diferencas finitas recuada 2 : " << rdiff2 << std::endl;
	std::cout <<  "Resultado por diferencas finitas centrada 1 : " << cdiff1 << std::endl;
	std::cout <<  "Resultado por diferencas finitas centrada 2 : " << cdiff2 << std::endl;
}

double f(double x){
	//return -0.1*(x*x*x*x) - 0.15*(x*x*x) - 0.5*(x*x) - 0.25*x + 1.2;
	//return std::exp(x);
	return x*std::log(x);
}
//Diferenças finitas avançada:
double fdiff_a_1st(double x, double h){
	double res = (f(x+h) - f(x)) / h;
	return res;
}
double fdiff_a_1st_v2(double x, double h){
	double res = (-f(x+2*h) + 4*f(x+h) - 3*f(x)) / (2*h);
	return res;
}
double fdiff_a_2nd(double x, double h){
	double res = (f(x+2*h) - 2*f(x+h) + f(x)) / (h*h);
	return res;
}
double fdiff_a_2nd_v2(double x, double h){
	double res = (-f(x+3*h) - 4*f(x+2*h) - 5*f(x+h) + 2*f(x)) / (h*h);
	return res;
}

//Diferenças finitas recuada:
double fdiff_r_1st(double x, double h){
	double res = (f(x) - f(x-h)) / h;
	return res;
}
double fdiff_r_1st_v2(double x, double h){
	double res = (3*f(x) - 4*f(x-h) + f(x-2*h)) / (2*h);
	return res;
}
double fdiff_r_2nd(double x, double h){
	double res = (f(x) - 2*f(x-h) + f(x-2*h)) / (h*h);
	return res;
}
double fdiff_r_2nd_v2(double x, double h){
	double res = (2*f(x) - 5*f(x-h) + 4*f(x-2*h) + f(x-3*h)) / (h*h);
	return res;
}
//Diferenças finitas centrada:
double fdiff_c_1st(double x, double h){
	double res = (f(x+h) - f(x-h)) / (2*h);
	return res;
}
double fdiff_c_1st_v2(double x, double h){
	double res = (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h)) / (12*h);
	return res;
}
double fdiff_c_2nd(double x, double h){
	double res = (f(x+h) - 2*f(x) + f(x-h)) / (h*h);
	return res;
}
double fdiff_c_2nd_v2(double x, double h){
	double res = (-f(x+2*h) + 16*f(x+h) - 30*f(x) + 16*f(x-h) - f(x-2*h)) / (12*(h*h));
	return res;
}