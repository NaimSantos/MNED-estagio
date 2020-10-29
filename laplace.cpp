#include <iostream>
#include <cmath>

#define P 4

typedef float newvar[P+1][P+1];
void inputrow(int i, newvar u){

	std::cout << "Informe o valor de u[" << i << ",j], j=1," << P;
    for(int j=1; j<=P; j++)
		std::cin >> u[i][j];

}
void inputcol(int j, newvar u){

	std::cout << "Informe o valor de u[i," << j << "], i=2," << P-1;
    for(int i=2; i<=P-1; i++)
		std::cin >> u[i][j];

}
void oput(newvar u, int wd, int prsn){

    for(int i=1; i<=P; i++){
		for(int j=1; j<=P; j++){
		std::cout << wd << ',' << prsn << ',' << u[i][j] << '\n';
		}
	}

}
int main(int argc, char* arg[]){
	newvar u;
	double mer, e, t;

	for(int i=1; i<=P; i++)
		for(int j=1; j<=P; j++)
			u[i][j]=0;
	printf("\nInforme a condicao de contorno:\n");
	inputrow(1, u);
	inputrow(P, u);
	inputcol(1, u);
	inputcol(P, u);
	
	const double ar = 0.001;	//erro
	const unsigned maxitr = 30;	//maximo de iterações
	for(int itr=1; itr<=maxitr; itr++){
			mer=0.0;
			for(int i=2; i<P-1; i++){
				for(int j=2; j<=P-1; j++){
					t = (u[i-1][j] + u[i+1][j] + u[i][j+1] + u[i][j-1]) / 4;
					e = std::abs(u[i][j] - t);
					if(e>mer)
						mer = e;
					u[i][j] = t;
				}
				std::cout << " Iteracao: " << itr << '\n';
				oput(u, 9, 2);
				if(mer<=ar){
					std::cout <<" Solucao:\n";
					oput(u, 8, 1);
					return 0;
				}
			}
		std::cout << "Numero de iteracoes insuficiente";
		return 1;
	}
}
