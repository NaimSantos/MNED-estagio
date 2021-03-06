#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

/*TO DO:
	usar vetores temporários c_temp e d_temp, para não escrever em c e d
	usar f como o vetor de saída.
*/
using std::vector;

void tdma_solver(const vector<double>& a, const vector<double>& b, vector<double>& c, vector<double>& d){
	auto n = static_cast<int>(d.size()-1);

	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++){
		c[i] = (c[i] ) / (b[i] - a[i]*c[i-1]);
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
	}

	d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

	for (int i = n; i-- > 0;){
		d[i] = d[i] - (c[i]*d[i+1]);
	}
}
int main(int argc, char* arg[]){
	vector<double> a1 = {0, 1, 2, 3};    //diagonal inferior
	vector<double> b1 = {2, 3, 5, 8};    //diagonal principal
	vector<double> c1 = {1, 2, 1, 0};    //diagonal superior
	
	vector<double> d1 = {7, 19, 31, 52};

	tdma_solver(a1, b1, c1, d1);

	std::cout << "T = ";
	for (auto &e: d1){
		std::cout << e << ' ';
	}
}
