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

  g(0) = 0; g(1) = 0; g(2) = 0; g(4) = 0; g(5) = 0; g(8) = 0;
}

void hessLJ(column_vector cluster, double rho, double E, int N, matrix<double>& H) {
  //compute the Hessian of particles under lennard-jones potential

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

              //compute term1
              double r1 = pow(r, -rho-2) * (1-pow(r,-rho));
              double T1 = (delta(i,j)-delta(i-DIMENSION*(particle1-particle2),j))*r1;
              
              //compute term 2
              double P2; double P3;
              if (particle1 == j / DIMENSION) {
                P2 = -(rho+2)*pow(r,-rho-4) * (cluster(j) - cluster(j - DIMENSION*(particle1-particle2)));
                P3 = rho * pow(r,-rho-2) * (cluster(j) - cluster(j - DIMENSION*(particle1-particle2)));
              }
              else{
                P2 = 0; P3 = 0;
                if (particle2 == j / DIMENSION) {
                  P2 = -(rho+2)*pow(r,-rho-4) * (cluster(j) - cluster(j - DIMENSION*(particle2-particle1)));
                  P3 = rho * pow(r,-rho-2) * (cluster(j) - cluster(j - DIMENSION*(particle2-particle1)));
                }
              }
              double T2 = (cluster(i)-cluster(i-DIMENSION*(particle1-particle2)))*P2*(1-pow(r,-rho));

              //compute term 3
              double T3 = (cluster(i)-cluster(i-DIMENSION*(particle1-particle2)))*P3*pow(r,-rho-2);
            
              //combine into a total
              S += T1 + T2 + T3;
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
