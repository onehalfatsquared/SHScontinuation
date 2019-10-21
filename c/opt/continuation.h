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

#define LOW 7.5
#define MED 7.5
#define HIGH 7.5

#define STEP 0.01


void testCV(double* clusters);
void descent(int N, int num_clusters, double* clusters, int sticky, int potential);
void parallelDescent(int N, int num_clusters, double* clusters, int sticky, int potential);


//functions to access shs and other clusters
int getNumClusters(int N); 
double getEnergy0(int sticky);
void getSHS(int N, int num_clusters, double* clusters);
void getFinite(int N, int sticky, int potential, double range, double* clusters);
void storeCluster(int N, column_vector out, int num, double* clusters);
void getCluster(int N, double* clusters, int num, column_vector& out);
void printClusters(int N, int num_clusters, double* clusters, double range, 
									 int sticky, int potential);
bool testSame(int N, double* clusters, int c1, int c2, double& distance);

//merges
bool isMerged(int cluster, std::vector<int> merged); 
void findMerges(int N, int num_clusters, int sticky, int potential); 
double rMin(int N, int cNum, double* clusters, double rho); 
void rMinFill(int N, int num_clusters, double* clusters, double* rMins, double rho);

