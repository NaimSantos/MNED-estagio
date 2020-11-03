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

void thomas_solver(const std::vector<double>& a, const std::vector<double>& b,
	const std::vector<double>& c,  const std::vector<double>& d, std::vector<double>& f)
{
	size_t N = d.size();


	std::vector<double> c_start(N, 0.0);
	std::vector<double> d_start(N, 0.0);
 
	c_start[0] = c[0] / b[0];
	d_start[0] = d[0] / b[0];

 
	for (int i=1; i<N; i++) {
		double m = 1.0 / (b[i] - a[i] * c_start[i-1]);
		c_start[i] = c[i] * m;
		d_start[i] = (d[i] - a[i] * d_start[i-1]) * m;
	}

                                                                                                                                              
	for (int i=N-1; i-- > 0; ) {
		f[i] = d_start[i] - c_start[i] * d[i+1];
	}
}


#endif //HEADER_H
