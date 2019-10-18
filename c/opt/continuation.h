#pragma once
#include "physics.h"
#include "dlib/optimization.h"

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


void testCV(double* clusters);



//functions to access shs and other clusters
int getNumClusters(int N); 
void getSHS(int N, int num_clusters, double* clusters);
void storeCluster(int N, column_vector out, int num, double* clusters);
void getCluster(int N, double* clusters, int num, column_vector& out);