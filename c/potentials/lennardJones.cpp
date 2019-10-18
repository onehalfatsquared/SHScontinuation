#include <math.h>
#include <algorithm>
#include "physics.h"

double ljP(double r, double rho, double E){
  //evaluate the lennard-jones potential derivative
  double Y = pow(1/r, rho);
  return  E * Y * (Y-2);
}

double ljPD(double r, double rho, double E){
  //evaluate the lennard-jones potential derivative
  double Y = pow(1/r, rho);
  return -2 * E * rho / r * Y * (Y-1);
}

double ljEval(double* particles, double rho, double E, int N) {
  //compute total energy of system with pairwise morse potential 
	double S = 0; 
  for (int i = 0; i < N; i++){
    for (int j = i+1; j < N; j++){
      if (j != i){
        double* Z = new double[DIMENSION];
        double r = euDist(particles, i, j, N, Z);
        S += ljP(r, rho, E); 
        delete []Z;
      }
    }
  }
  return S;
}

void ljGrad(double* particles, double rho, double E, int N, column_vector& g) {
  //compute gradient of energy of system with pairwise morse potential
	double* S = new double[DIMENSION];
	for (int i = 0; i < N; i++){
    S[0] = S[1] = 0;
#if (DIMENSION == 3)
    S[2] = 0;
#endif

    for (int j = 0; j < N; j++){
      if (j != i){
        double* Z = new double[DIMENSION];
        double r = euDist(particles, i, j, N, Z);
        for(int k = 0; k < DIMENSION; k++) {
          S[k] = S[k] + ljPD(r, rho, E) * Z[k];
        }
      	delete []Z;
    	}
		}
    for(int k = 0; k < DIMENSION; k++) {
      g(DIMENSION*i+k) = S[k];
    }
  }
  delete []S;
}

void hessLJ(double* particles, double rho, double E, int N, double* H) {
  //compute the Hessian of particles under morse potential


}
