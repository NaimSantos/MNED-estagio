#ifndef GENERIC_ALG_H
#define GENERIC_ALG_H

#include <cmath>

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

//print a bidimensional array of elements of type T
template <typename T>
void printarray2D(T** M, const int nrow, const int ncol){
	for(int i=0; i<nrow; i++){
		for(int j=0; j<ncol; j++){
			std::cout << M[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

#endif //HEADER_H