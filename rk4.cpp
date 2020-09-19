#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

double dnx(double, double [], double [], int );
double rk4_dn1(double(*)(double, double [], double [], int ), double, double, double [], double [], int);

int main(){
	const int n=4;                   // number of first-order equations 
	double ti, tf, dt, tmax;
	double xi[n], xf[n];
	double v0, a0;
	int i, key;
	static const double g = 9.81;               // free fall acceleration in m/s^2
	static const double m = 1.0;                // mass of a projectile in kg
	static const double rad = 3.1415926/180.0;  // radians

    std::ofstream file;
    file.open ("table01c.dat");

    file.precision(3);
    file.setf(ios::scientific | ios::showpoint);
    std::cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);

    ti = 0.0;                // initial value for variable t
    v0 = 180.0;              // initial speed (m/s)
    a0 =  45.0;              // initial angle (degrees)
    xi[0] = 0.0;             // initial position in x (m)
    xi[1] = 0.0;             // initial position in y (m)
    xi[2] = v0*cos(a0*rad);  // initial speed in x direction (m.s)
    xi[3] = v0*sin(a0*rad);  // initial speed in y direction (m/s)
    
    dt = 0.2;             // step size for integration (s)
    tmax = 60.0;          // integrate till tmax (s)
	
    file << std::setw(12) << "t"  << std::setw(12) << "x"   << std::setw(12) << "y"
         << std::setw(12) << "x'" << std::setw(12) << "y'"  << endl; 

/* integration of ODE */
    while (ti <= tmax){

     file << std::setw(12) << ti   << std::setw(12) << xi[0] << std::setw(12) << xi[1]
          << std::setw(12) << xi[2]<< std::setw(12) << xi[3] << endl;

     if (xi[1] < 0.0) break; 

     tf = ti + dt;
     rk4_dn1(dnx, ti, tf, xi, xf, n);

        ti = tf;
        for (i = 0; i<=n-1; i = i+1){ 
           xi[i] = xf[i];
        }
   }
    system ("pause");
    return 0;
}

/*============================================================== 
  System of first order differential equations for the RK solver
  
  For a system of n first-order ODEs
  x [] array - x values
  dx[] array - dx/dt values
  
  For a system of n/2 second order ODEs follow the agreement
  In:  x[] array 
  # first n/2 elements are x
  # last  n/2 elements are dx/dt
  Out: dx[] array
  # first n/2 elements are first order derivatives  (or dx/dt)
  # last  n/2 elements are second order derivatives (or d2x/dt2)
  example: 2D projectile motion in (x,y) plane
  In           Out
  x[0] - x     dx[0] - x'
  x[1] - y     dx[1] - y'
  x[2] - x'    dx[2] - x"
  x[3] - y'    dx[3] - y"
==============================================================*/

double dnx(double t, double x[], double dx[], int n){
	/* first order */
		dx[0] = x[2];
		dx[1] = x[3];
	/* second order */
		dx[2] = 0.0;
		dx[3] = (-1.0)*g;
}    
double rk4_dn1(double(dnx)(double, double [], double [], int),double ti, double tf, double xi[], double xf[], int n){
      double h, t, x[n], dx[n];
      double k1[n],k2[n],k3[n],k4[n];
      int j;
      h = tf-ti;
      t = ti;
      dnx(t, xi, dx, n);
      for (j = 0; j<=n-1; j = j+1) {
          k1[j] = h*dx[j];
          x[j]  = xi[j] + k1[j]/2.0;  
        }      
      dnx(t+h/2.0, x, dx, n);
      for (j = 0; j<=n-1; j = j+1){
          k2[j] = h*dx[j];
          x[j]  = xi[j] + k2[j]/2.0;  
        }
      dnx(t+h/2.0, x, dx, n);
      for (j = 0; j<=n-1; j = j+1){
          k3[j] = h*dx[j];
          x[j]  = xi[j] + k3[j];  
        }   
      dnx(t+h, x, dx, n);
      for (j = 0; j<=n-1; j = j+1){
          k4[j] = h*dx[j];
          xf[j] = xi[j] + k1[j]/6.0+k2[j]/3.0+k3[j]/3.0+k4[j]/6.0;
        }      
    return 0.0;
}