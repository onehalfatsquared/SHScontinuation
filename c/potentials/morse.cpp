#include <math.h>
#include <algorithm>
#include "physics.h"

double morseP(double r, double rho, double E) {
  //evaluate the morse potential 
  double Y = exp(-(rho) * ((r) - 1));
  return E * ( Y*Y - 2*Y);
}

double morsePD(double r, double rho, double E) {
  //evaluate the morse potential derivative
  double Y = exp(-(rho) * ((r) - 1));
  return -2 * (rho) * (E) * ( Y*Y -Y);
}


double morseEval(double* particles, double rho, double E, int N) {
  //compute total energy of system with pairwise morse potential 
	double S = 0; 
  for (int i = 0; i < N; i++) {
    for (int j = i+1; j < N; j++) {
      if (j != i){
        double* Z = new double[DIMENSION];
        double r = euDist(particles, i, j, N, Z);
        S += morseP(r, rho, E); 
        delete []Z;
      }
    }
  }
  return S;
}

void morseGrad(double* particles, double rho, double E, int N, column_vector& g) {
  //compute gradient of energy of system with pairwise morse potential
	double* S = new double[DIMENSION];
	for (int i = 0; i < N; i++) {
    S[0] = S[1] = 0;
#if (DIMENSION == 3)
    S[2] = 0;
#endif

  for (int j = 0; j < N; j++) {
    if (j != i) {
      double* Z = new double[DIMENSION];
      double r = euDist(particles, i, j, N, Z);
      for (int k = 0; k < DIMENSION; k++) {
        S[k] = S[k] + morsePD(r, rho, E) * Z[k];
      }
    	delete []Z;
  	}
	}
  for (int k = 0; k < DIMENSION; k++)
    g(DIMENSION*i+k) = S[k];
  }
  delete []S;

  g(0) = 0; g(1) = 0; g(2) = 0; g(4) = 0; g(5) = 0; g(8) = 0;
}

void hessMorse(double* particles, double rho, double E, int N, double* H) {
  //compute the Hessian of particles under morse potential


}
