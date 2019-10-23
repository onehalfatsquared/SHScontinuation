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

bool isDoF(int i) {
  //return true if i is a true degree of freedom
  // return false if i is 0,1,2,4,5,8

  if (i == 0 || i == 1 || i == 2 || i == 4 || i == 5 || i == 8) {
    return false;
  }

  return true;
}

int delta(int i, int j) {
  //kronecker delta fn
  if (i == j) {
    return 1;
  }
  return 0;
}

void hessMorse(column_vector cluster, double rho, double E, int N, matrix<double>& H) {
  //compute the Hessian of particles under morse potential

  int entries = DIMENSION*N;

  double* particles = new double[entries];
  c2p(cluster, particles, N);

  double* Z = new double[DIMENSION];

  for (int i = 0; i < entries; i ++) {
    if (isDoF(i)) {
      for (int j = 0; j < entries; j++) {
        if (isDoF(j)) { //if we get here, compute contribution from this entry
          //init the sum
          double S = 0;

          //get the particle of the i index
          int particle1 = i / DIMENSION; 

          //loop over all other particles
          for (int particle2 = 0; particle2 < N; particle2++) {
            if (particle2 != particle1) {
              //get quantities to use in later expressions
              double r = euDist(particles, particle1, particle2, N, Z);
              double Y = exp(-rho*(r-1));

              //compute term1
              double T1 = ((delta(i,j)-delta(i-DIMENSION*(particle1-particle2),j))/r) * Y*(1-Y);
              
              //compute term 2
              double P2;
              if (particle1 == j / DIMENSION) {
                P2 = -1.0/(r*r*r) * (cluster(j) - cluster(j - DIMENSION*(particle1-particle2)));
              }
              else{
                P2 = 0;
                if (particle2 == j / DIMENSION) {
                  P2 = -1.0/(r*r*r) * (cluster(j) - cluster(j - DIMENSION*(particle2-particle1)));
                }
              }
              double T2 = (cluster(i)-cluster(i-DIMENSION*(particle1-particle2)))*Y*(1.0-Y)*P2;

              //compute term 3
              double P3 = rho * Y * r*r * P2;
              double T3 = ((cluster(i)-cluster(i-DIMENSION*(particle1-particle2))) / r) * (1.0-Y) * P3;

              //compute term 4
              double T4 = ((cluster(i)-cluster(i-DIMENSION*(particle1-particle2))) / r) * Y * (-P3);
            
              //combine into a total
              S += T1 + T2 + T3 + T4;
            }
          }
          //add entry into matrix 
          H(i,j) = 2 * E * rho * S;
        }
      }
    }
  }

  delete []Z; delete []particles;


}
