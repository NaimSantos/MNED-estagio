#ifndef GENERIC_ALG_H
#define GENERIC_ALG_H

#include <cmath>
#include <vector>
#include <iostream>


//Uma estimativa para o valor de PI:
constexpr double NPI = 4 * std::atan(1);


//Conversão de graus para radianos:
double torad(double value){
	return value * NPI / 180 ;
}


//Retorna se um arranjo bidimensional é uma matriz diagonalmente dominante:
bool isdiagonaldom(double** M, const int nrows, const int ncol){

	for (int i = 0; i < nrows; i++){

		double sum = 0.0;
		double diag = std::abs(M[i][i]);

		for (int j = 0; j < ncol; j++){
			if (i==j)
				continue;
			sum += std::abs(M[i][j]);
		}
		if (sum > diag)
			return false;
	}
	return true;
}


//Printa um arranjo bidimensional de elementos de um dado tipo T:
template <typename T>
void printarray2D(T** M, const int nrow, const int ncol){
	for(int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			std::cout << M[i][j] << ' ';
		}
		std::cout << '\n';
	}
}


//TDMA: TriDiagonal Matrix Algorithm ("Thomas")
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b,
	const std::vector<double>& c,  const std::vector<double>& d, std::vector<double>& f)
{
	int N = d.size();

	std::vector<double> c_start(N, 0.0);
	std::vector<double> d_start(N, 0.0);

	c_start[0] = c[0] / b[0];
	d_start[0] = d[0] / b[0];

	for (int i=1; i<N; i++){
		double m = 1.0 / (b[i] - a[i] * c_start[i-1]);
		c_start[i] = c[i] * m;
		d_start[i] = (d[i] - a[i] * d_start[i-1]) * m;
	}

	for (int i=N-1; i-- > 0; ) {
		f[i] = d_start[i] - c_start[i] * d[i+1];
	}
}


#endif //HEADER_H
