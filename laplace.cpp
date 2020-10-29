#include <iostream>
#include <cmath>
#include "alg_utilities.h"

void inputrow(int i, double** T, const int dim);
void inputcol(int j, double** T, const int dim);
void oput(double** T, int wd, int prsn, const int dim);

int main(int argc, char* arg[]){

	//Aloca a matriz da malha:
	const int dim = 3;
	double** T{new double*[dim] {nullptr}};
	for (int i = 0; i < dim; ++i)
		T[i] = new double[dim] {0.0};
	
	printarray2D<double>(T, dim, dim);


	std::cout << "\nInforme a condicao de contorno:\n";
	inputrow(1, T, dim);
	inputrow(dim, T, dim);
	inputcol(1, T, dim);
	inputcol(dim, T, dim);

	const double ar = 0.001;	//erro
	const unsigned maxitr = 30;	//maximo de iterações
	double temp = 0.0;
	double mer = 0.0;
	double erro = 0.0;
	
	for(int itr=1; itr<=maxitr; itr++){
		for(int i=2; i<dim-1; i++){
			for(int j=2; j<=dim-1; j++){
				temp = (T[i-1][j] + T[i+1][j] + T[i][j+1] + T[i][j-1]) / 4;
				erro = std::abs(T[i][j] - temp);
				if(erro>mer)
					mer = erro;
				T[i][j] = temp;
			}
			std::cout << " Iteracao: " << itr << '\n';
			oput(T, 9, 2, dim);
			if(mer<=ar){
				std::cout <<" Solucao:\n";
				oput(T, 8, 1, dim);
			}
		}
	}
	//Desalocar a matrix:
	for (int i = 0; i < dim; ++i)
		delete[] T[i];
	return 0;
}
void inputrow(int i, double** T, const int dim){

	std::cout << "Informe o valor de T[" << i << ",j], j=1," << dim << ": ";
	for(int j=1; j < dim; j++)
		std::cin >> T[i][j];

}
void inputcol(int j, double** T, const int dim){

	std::cout << "Informe o valor de T[i," << j << "], i=2," << dim-1 << ": ";
	for(int i=2; i < dim-1; i++)
		std::cin >> T[i][j];

}
void oput(double** T, int wd, int prsn, const int dim){

	for(int i=1; i < dim; i++){
		for(int j=1; j<dim; j++)
			std::cout << wd << ',' << prsn << ',' << T[i][j] << '\n';
	}

}