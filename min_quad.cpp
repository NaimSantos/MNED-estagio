#include <iostream>
#include <cmath>
#include <array>
#include <iomanip>

struct LinReg
{
	double a;
	double b;
	double r;
	double eps;
};

template <typename T>
T average(T* A, const unsigned int n);

template <typename T>
double desvio(T* A, const unsigned int n, const T aver);

template <typename T1, typename T2>
double covar(T1* A1, T2* A2, const unsigned int n);

template<typename T1, typename T2>
LinReg RegressaoLinear(T1* A1, T2* A2, const unsigned int n);

int main(int argc, char* arg[])
{
	int X[7] {1, 2, 3, 4, 5, 6, 7};
	double Y[7] {2, 4.1, 5.8, 8.3, 10, 12.3, 14.3};
	auto resultado = RegressaoLinear<int, double>(X, Y, std::size(X));
	
	std::cout << std::setprecision(10) << "a = " << resultado.a << "\b = " << resultado.b
		<< "\r = " << resultado.r << "\eps = " << resultado.eps << std::endl;
}

template <typename T>
T average(T* A, const unsigned int n){
	T sum {};
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
//Retorna um struct com 4 membros, .a, .b, .r, .eps, coeficientes 
//.a= inclinação, .b= intersseção, .r=coeficiente de Pearson e .eps=erro padrao de y;
template<typename T1, typename T2>
LinReg RegressaoLinear(T1* A1, T2* A2, const unsigned int n){

	LinReg lin_reg;

	auto covarxy = covar<T1, T2>(A1, A2, n);
	auto xmed = average<T1>(A1, n);
	auto ymed = average<T2>(A2, n);
	auto desvx = desvio<T1>(A1, n, xmed);
	auto desvy = desvio<T2>(A2, n, ymed);

	//a
	auto a = (covarxy / (desvx*desvx));
	lin_reg.a = a;

	//b
	auto b = ymed - (a * xmed);
	lin_reg.b = b;

	//r
	auto r = covarxy /( desvx * desvy);
	lin_reg.r = r;

	//epsilon_y
	{
		double sum = 0.0;
		for (int i = 0; i < n; i++){
			sum = std::pow((A2[i] - ((a*A1[i]) + b)), 2) / (n-2);
		}
		auto eps = std::sqrt(sum);
		lin_reg.eps = eps;
	}

	return lin_reg;
}
