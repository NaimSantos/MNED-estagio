#include <iostream>
#include <iomanip>
#include <cmath>

void GS_Solver(double** A, const int m, const int n, double* B, const unsigned int n_eq, const float eps, double* X);

int main(){

	//Alocação: A[dim1][dim2], B[dim1], C[dim1];
	const int dim1 = 3;
	const int dim2 = 3;

	double** A{new double*[dim1] {}};
	for (int i = 0; i < dim1; ++i)
		A[i] = new double[dim2] {};

	double* B = new double [dim1];
	double* X = new double [dim1];

	//Exemplo
	A[0][0] = 4.0; A[0][1] = 1.0; A[0][2] = -1.0;
	A[1][0] = 2.0; A[1][1] = 7.0; A[1][2] = 1.0;
	A[2][0] = 1.0; A[2][1] = -3.0; A[2][2] = 12.0;

	B[0]= 3.0; B[1]=19.0; B[2] = 31.0;
	X[0]= 0.0; X[1]=0.0; X[2]=0.0;

	const double eps = 0.00001; //define o erro permitido, usado como criterio de parada

	GS_Solver(A, dim1, dim2, B, dim1, eps, X);

	//Desalocar as matrizes:
	for (int i = 0; i < dim1; ++i)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;

	return 0;
}

void GS_Solver(double** A, const int m, const int n, double* B, const unsigned int n_eq, const float eps, double* X){

	double* Y = new double [n_eq];	//matrix auxiliar
	double* E = new double [n_eq]; //necessária para estimar o erro de uma iteração a outra
	for (int i = 0; i < n_eq; i++)
		E[i] = X[i];
	
	int counter = 1; //Contar iterações apenas pro caso da tolerancia nao ser atingida. Se a matriz é diagonal dominante, a convergência é garantida
	bool teste = false;

	while(!teste && counter<50){
		std::cout << "Iteracao " << std::setprecision(10) << counter << std::endl;
		for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j==i){
					continue;
				}
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				X[i] = Y[i]; //escreve em X a estimativa encontrada
				
			}
			auto res = std::abs(((X[i] - E[i]) / X[i])) <= eps;
			(!res) ? teste = false : teste = true;
			std::cout<< "x" << i + 1 << " = " << Y[i] << std::endl;
			E[i]=X[i];
		}
		counter++;
		std::cout << std::endl;
	}
	
	delete[] E;
	delete[] Y;
}
