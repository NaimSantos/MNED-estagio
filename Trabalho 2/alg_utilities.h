#ifndef GENERIC_ALG_H
#define GENERIC_ALG_H

#include <cmath>
#include <vector>
#include <iostream>

constexpr double NPI = 4 * std::atan(1);
//coverts degree to radians:
double torad(double value){
	return value * NPI / 180 ;
}
//returns if a bidimensional array is diagonally domminant matrix
bool isdiagonaldom(double** M, const int nrows, const int ncol){

	for (int i = 0; i < nrows; i++){

		double sum = 0.0;
		double dia = std::abs(M[i][i]);

		for (int j = 0; j < ncol; j++){
			if (i==j)
				continue;
			sum += std::abs(M[i][j]);
		}
		if (sum > dia)
			return false;
	}
	return true;
}

//prints a bidimensional array of elements of type T
template <typename T>
void printarray2D(T** M, const int nrow, const int ncol){
	for(int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			std::cout << M[i][j] << ' ';
		}
		std::cout << '\n';
	}
}
/*
//TDMA: TriDiagonal Matrix Algorithm ("Thomas")
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, const std::vector<double>& d,
	std::vector<double>& f, std::vector<double>& c_temp, std::vector<double>& d_temp){
	auto N = d.size();

	std::fill(c_temp.begin(), c_temp.end(), 0.0);
	std::fill(d_temp.begin(), d_temp.end(), 0.0);

	c_temp[0] = c[0] / b[0];
	d_temp[0] = d[0] / b[0];

	for (int i = 1; i < N; i++){
		double m = 1.0 / (b[i] - a[i] * c_temp[i-1]);
		c_temp[i] = c[i] * m;
		d_temp[i] = (d[i] - a[i] * d_temp[i-1]) * m;
	}
	for (int i = N-1; i > 0; i--){
		f[i] = d_temp[i] - c_temp[i] * d[i+1];
	}
}

*/
#endif //HEADER_H
