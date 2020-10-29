#ifndef GENERIC_ALG_H
#define GENERIC_ALG_H

#include <cmath>

//returns if a bidimensional array is diagonal domminant matrix
bool isdiagonaldom(double** M, const int nrows, const int ncol){

	for (int i = 0; i < nrows; i++){

		double sum = 0.0;
		double dia = std::abs(M[i][i]);

		for (int j = 0; j < ncol; j++){
			if (i==j)
				continue;
			sum += std::abs(M[i][j]);
		}
		if (sum >= dia)
			return false;
	}
	return true;
}

#endif //HEADER_H