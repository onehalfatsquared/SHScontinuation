#pragma once
#include "physics.h"
#include "dlib/optimization.h"
#include <omp.h>
#include <vector>

class Parameters {
	public:
		Parameters(int N_, int potential_, double rho_, double E_);
		~Parameters();

		double getU(column_vector cluster) const;
		column_vector getGrad(column_vector cluster) const;

	private:
		int N, potential;
		double E, rho;

};

#define LOW 7.51
#define MED 8.3
#define HIGH 9.5

#define STEP 0.01

//test functions
void testCV(double* clusters);

//continuation/descent functions
double getEnergy0(int sticky);
void descent(int N, int num_clusters, double* clusters, int sticky, int potential);
void parallelDescent(int N, int num_clusters, double* clusters, int sticky, int potential);
int countNegs(column_vector eig, int N);
bool isMin(eigenvalue_decomposition<matrix<double>> Hd, int N, std::vector<int>& indices);
void reOpt(column_vector& X, double range, double E, int N, int potential, std::vector<int> indices,
					 Parameters p, eigenvalue_decomposition<matrix<double>> Hd); 


//functions to access shs and other clusters
int getNumClusters(int N); 
void getSHS(int N, int num_clusters, double* clusters);
void getFinite(int N, int sticky, int potential, double range, double* clusters);
void storeCluster(int N, column_vector out, int num, double* clusters);
void getCluster(int N, double* clusters, int num, column_vector& out);
void printClusters(int N, int num_clusters, double* clusters, double range, 
									 int sticky, int potential);

//functions to detect merges and compute properties of clusters
bool isMerged(int cluster, std::vector<int> merged); 
void findMerges(int N, int num_clusters, int sticky, int potential); 
double rMin(int N, int cNum, double* clusters, double rho); 
void rMinFill(int N, int num_clusters, double* clusters, double* rMins, double rho);
bool testSame(int N, double* clusters, int c1, int c2, double& distance); 
bool testSame(int N, column_vector c1, column_vector c2, double& distance);
double computeClusterMetric(int N, column_vector c1, column_vector c2, int which);

