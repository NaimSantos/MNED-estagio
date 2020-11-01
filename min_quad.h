#ifndef LIN_REG_H
#define LIN_REG_H

#include <cmath>
#include <array>

struct LinReg
{
	double a;
	double b;
	double r;
	double eps;
};

template <typename T>
double average(T* A, const unsigned int n);

template <typename T>
double desvio(T* A, const unsigned int n, const T aver);

template <typename T1, typename T2>
double covar(T1* A1, T2* A2, const unsigned int n);

template<typename T1, typename T2>
LinReg RegressaoLinear(T1* A1, T2* A2, const unsigned int n);

template <typename T>
double average(T* A, const unsigned int n){
	double sum {};
	for (int i = 0; i < n ; i++)
		sum += A[i];
	return sum / n;
}

template <typename T>
double desvio(T* A, const unsigned int n, const T aver){
	double sum = 0.0;
	for (int i = 0; i < n ; i++){
		sum += std::pow(A[i] - aver, 2) / (n-1);
	}
	return std::sqrt(sum);
}

template <typename T1, typename T2>
double covar(T1* A1, T2* A2, const unsigned int n){
	auto x_med = average<T1>(A1, n);
	auto y_med = average<T2>(A2, n);
	double sum = 0.0;

	for (int i = 0; i < n ; i++){
		sum += ((A1[i] - x_med)*(A2[i] - y_med)) / (n-1);
	}
	return sum;
}
//Retorna um struct com 4 componentes, coeficientes da regressão: a, b, r e eps
template<typename T1, typename T2>
LinReg RegressaoLinear(T1* A1, T2* A2, const unsigned int n){

	LinReg res;

	auto covarxy = covar<T1, T2>(A1, A2, n);
	auto xmed = average<T1>(A1, n);
	auto ymed = average<T2>(A2, n);
	auto desvx = desvio<T1>(A1, n, xmed);
	auto desvy = desvio<T2>(A2, n, ymed);

	//Inclinação:
	auto a = (covarxy / (desvx*desvx));
	res.a = a;

	//Interseção:
	auto b = ymed - (a * xmed);
	res.b = b;

	//coeficiente de correlação de Pearson
	auto r = covarxy /( desvx * desvy);
	res.r = r;

	//erro padrão em y:
	double sum = 0.0;
	for (int i = 0; i < n; i++){
		sum += std::pow((A2[i] - ((a*A1[i]) + b)), 2) / (n-2);
	}
	auto eps = std::sqrt(sum);
	res.eps = eps;


	return res;
}
#endif //LIN_REG_H
