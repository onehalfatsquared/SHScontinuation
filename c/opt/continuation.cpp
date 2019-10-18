#include "physics.h"
#include "continuation.h"
#include "dlib/optimization.h"
#include <iostream>
#include <fstream>
#include "string.h"

using namespace dlib;

Parameters::Parameters(int N_, int potential_, double rho_, double E_) {
	N = N_; potential = potential_; rho = rho_; E = E_;
}

Parameters::~Parameters() {};

double Parameters::getU(column_vector cluster) const {
	//return the energy for given set of parameters

	double particles[DIMENSION*N]; c2p(cluster, particles, N);
	double U;

	if (potential == 0) {
		U = morseEval(particles, rho, E, N);
	}
	else if (potential == 1) {
		U = ljEval(particles, rho, E, N);
	}

	return U;

}

column_vector Parameters::getGrad(column_vector cluster) const {
	//return the gradient

	double particles[DIMENSION*N]; c2p(cluster, particles, N);
	column_vector g(DIMENSION*N); 

	if (potential == 0) {
		morseGrad(particles, rho, E, N, g);
	}
	else if (potential == 1) {
		ljGrad(particles, rho, E, N, g);
	}

	return g;

}



void testCV(double* clusters) {
	column_vector X(18); getCluster(6,clusters, 1, X);
	column_vector Y(18); getCluster(6,clusters, 0, Y);

	Parameters p = Parameters(6, 0, 30, 1);

	//std::cout << derivative([&p](const column_vector& a) {return p.getU(a);},1e-8)(X); 
	//std::cout << p.getGrad(X);
	//std::cout << (p.getU(Y)-p.getU(X))/1e-4 << "\n";
	//std::cout << (p.getU(Y)) << "\n";
	std::cout << X << "\n";
	std::cout << Y << "\n";

	find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
             objective_delta_stop_strategy(1e-13).be_verbose(), // Stop when the change in rosen() is less than 1e-7
             [&p](const column_vector& a) {return p.getU(a);}, 
             [&p](const column_vector& a) {return p.getGrad(a);}, 
             X, -10000);
	find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
             objective_delta_stop_strategy(1e-13).be_verbose(), // Stop when the change in rosen() is less than 1e-7
             [&p](const column_vector& a) {return p.getU(a);}, 
             [&p](const column_vector& a) {return p.getGrad(a);}, 
             Y, -10000);
    // Once the function ends the starting_point vector will contain the optimum point 
    // of (1,1).
    std::cout << "rosen solution:\n" << X << "\n";
    std::cout << "rosen solution:\n" << Y << "\n";
}

int getNumClusters(int N) {
	//return the number of clusters for given N

	int tot = 0; 
	
	if (N == 6) tot = 2;
	else if (N == 7) tot =  5;
	else if (N == 8) tot = 13;
	else if (N == 9) tot = 52; 

	return tot;
}

void getCluster(int N, double* clusters, int num, column_vector& out) {
	//access the num cluster in array and store in col vec

	//check if this is actually a cluster
	int elements = DIMENSION*N;
	int total = getNumClusters(N);

	if (num > total-1) {
		fprintf(stderr, "Trying to access non-existant cluster %d of %d\n", num+1, total);
		return;
	}

	//store cluster in column vec
	for (int i = 0; i < elements; i++) {
		out(i) = clusters[num*elements+i];
	}
}

void storeCluster(int N, column_vector out, int num, double* clusters) {
	//store the values in the vector back in the array

	//check if this is actually a cluster
	int elements = DIMENSION*N;
	int total = getNumClusters(N);

	if (num > total-1) {
		fprintf(stderr, "Trying to access non-existant cluster %d of %d\n", num+1, total);
		return;
	}

	//store cluster in column vec
	for (int i = 0; i < elements; i++) {
		clusters[num*elements+i] = out(i);
	}
}

void getSHS(int N, int num_clusters, double* clusters) {
	//import the starting shs clusters from file

	//file location
	std::string filename = "input/n" + std::to_string(N) + "adjust.txt";

	//open the file
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return;
	}

	//store the entries
	double val; int entry = 0;
	while (in_str >> val) {
		clusters[entry] = val;
		std::cout << clusters[entry] << "\n";
		entry++;
	}
}
