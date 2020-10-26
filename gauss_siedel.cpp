#include <iostream> 
#include <cmath>
#include <fstream>

void GS_Solver(double** A, const int m, const int n, double* B, const unsigned int n_eq, const float eps, double* X);

int main(){ 

	const int dim1 = 2;
	const int dim2 = 2;

	double** A{new double*[dim1] {}};
	for (int i = 0; i < dim1; ++i)
		A[i] = new double[dim2] {};

	double* B = new double [dim1];
	double* X = new double [dim1];


	const double eps = 0.0001; //define o erro permitido

	//Preencher aqui as matrizes A e B:
	/*exemplos:
	for (int i = 0; i < dim1; i++ ){
		for (int j=0; j<dim2; j++){
			A[i][j]=i+j;
		}
	}
	for (int i = 0; i < dim1; i++ ){
		B[i]=i+1;
	}
	*/
	/*Exemplo de sistema:
		2x1+ x2 = 1
		3x1+ 4x2 = -1
	tolerancia: 10^-4
	vetor inicial: x= [0, 0]
	*/
	A[0][0] = 2; A[0][1] = 1;
	A[1][0] = 3; A[1][1] = 4;

	B[0]= 1.0; B[1]=-1.0;
	X[0]= 0.0; X[1]=0.0;

	GS_Solver(A, dim1, dim2, B, dim1, eps, X);

	//desalocar as matrizes:
	for (int i = 0; i < dim1; ++i)
		delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] X;

	return 0;

}
/*Parametros:
	A: matriz de coeficientes, no tamanho m x n
	m,n: numero de linhas e colunas de A, respectivamente
	B: matriz de termos independentes
	n_eq = numero de linhas de B
	eps: tolerancia
	X: matriz com as estimativas iniciais, onde serão escritos os resultados
*/
void GS_Solver(double** A, const int m, const int n, double* B, const unsigned int n_eq, const float eps, double* X){

	double Y[n_eq] = {};	//matrix auxiliar, necessária para estimar o erro de uma iteracao a outra

	int counter = 1;
	/*
		Contar iterações apenas pro caso da tolerancia nao ser atingida.
		No entanto, se a matriz é diagonal dominante, a convergência é garantida
	*/
	bool teste = false;

	do{
		std::cout << "Iteracao " << counter << std::endl;
		for (int i=0; i<m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j=0; j<n; j++){
				if (j==i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				teste = ((Y[i]-X[i])/Y[i]) < eps;
				X[i] = Y[i]; //escreve em X a estimativa encontrada
			}
			std::cout<< "x" << i + 1 << " = "<< Y[i] << std::endl;
		}
		counter++;
		std::cout << std::endl;
	}
	while(!teste && counter<20); 
}
