#include "physics.h"
#include "continuation.h"
#include "dlib/optimization.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include "string.h"

using namespace dlib;

/* Define the class object to store all the parameters for an optimization call */

/********************************************************************************/

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

/*************************************************************************/

/* Functions to read, write and edit the input and output clusters arrays */

/*************************************************************************/

int getNumClusters(int N) {
	//return the number of clusters for given N

	int tot = 0; 
	
	if (N == 6) tot = 2;
	else if (N == 7) tot =  5;
	else if (N == 8) tot = 13;
	else if (N == 9) tot = 52; 

	return tot;
}

double getEnergy0(int sticky) {
	double k;

	if (sticky == 0) k = LOW;
	else if (sticky == 1) k = MED;
	else if (sticky == 2) k = HIGH;

	return k;
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

void printClusters(int N, int num_clusters, double* clusters, double range, 
									 int sticky, int potential) {
	//print all the new clusters out to a file

	//file location
	std::string filename = "output/n" + std::to_string(N);
	if (potential == 0) filename += "rho";
	else if (potential == 1) filename += "m";
	std::string s = std::to_string(range); s.erase(s.size()-4,s.size());
	filename += s;
	if (sticky == 0) filename += "LOW.txt";
	else if (sticky == 1) filename += "MED.txt";
	else if (sticky == 2) filename += "HIGH.txt";

	std::ofstream ofile;
	ofile.open(filename);

	for (int val = 0; val < num_clusters*N*DIMENSION; val++) {
		ofile << std::fixed << std::setprecision(16) << clusters[val] << ' ';
	}
	
}

bool testSame(int N, double* clusters, int c1, int c2, double& distance) {
	//check if two clusters are the same just by comparing list of particle distances

	//set tolerance for sorted distances
	double tol = 1e-6;

	//store cluster sin column vectors
	column_vector cluster1(DIMENSION*N), cluster2(DIMENSION*N);
	getCluster(N, clusters, c1, cluster1); getCluster(N, clusters, c2, cluster2); 

	//make particle arrays
	double* particles1 = new double[N*DIMENSION];
	double* particles2 = new double[N*DIMENSION];
	c2p(cluster1, particles1, N); c2p(cluster2, particles2, N);

	//make vectors with distances
	std::vector<double> d1, d2;
	double* Z = new double[DIMENSION];
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			d1.push_back(euDist(particles1, i, j, N, Z)); 
			d2.push_back(euDist(particles2, i, j, N, Z)); 
		}
	}

	//sort the distance vectors
	sort(d1.begin(), d1.end()); sort(d2.begin(), d2.end()); 

	//compute summed distance between elements
	distance = 0;
	for (int i = 0; i < d1.size(); i++) {
		distance += abs(d1[i]-d2[i]);
		//std::cout << d1[i] << ' ' << d2[i] << ' ' << distance << "\n";
	}

	//free memory
	delete []particles1; delete []particles2; delete []Z;

	//return if they are the same
	if (distance < tol) return true;
	return false;

}


/****************************************************************************/

/* Functions to set up and do the optimizations */

/****************************************************************************/

void descent(int N, int num_clusters, double* clusters, int sticky, int potential) {
	//perform the descent with given parameters 

	double range = 50;               //set the initial range 
	double E0 = getEnergy0(sticky);  //energy at range 50
	double E = E0;                   //energy at arbitrary range
	double kappa, kappaD;            //kappa and its derivative
	stickyF(E0, range, 1, 0, kappa, kappaD); //get initial kappa
	double end = 1;

	range = range - STEP;
	while (range > end) {
		//get the energy from sticky parameter value 
		E = stickyNewton(E0, range, kappa, 1); 

		//construct the parameters object 
		Parameters p = Parameters(N, potential, range, E);

		//loop over the clusters 
		for (int c = 0; c < num_clusters; c++) {
			//get the previous cluster
			column_vector X(DIMENSION*N); 
			getCluster(N, clusters, c, X);

			//perform minimization
			find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
             objective_delta_stop_strategy(1e-13), //.be_verbose(), // Stop when the change in rosen() is less than 1e-7
             [&p](const column_vector& a) {return p.getU(a);}, 
             [&p](const column_vector& a) {return p.getGrad(a);}, 
             X, -10000);

			//store the new minimum
			storeCluster(N, X, c, clusters);
		}

		//output the cluters to a file
		printClusters(N, num_clusters, clusters, range, sticky, potential);

		//test if same 
		double d;
		std::cout << range << ' ' << testSame(N, clusters, 0, 1, d) << ' ' << d << "\n";

		//reduce the range parameter
		range -= STEP;
		
	}

	
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


