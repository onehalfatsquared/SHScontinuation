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
		//std::cout << clusters[entry] << "\n";
		entry++;
	}
}

void getFinite(int N, int sticky, int potential, double range, 
							 double* clusters) {
	//import the clusters found at given range from descent

	//file location
	std::string filename = "output/n" + std::to_string(N);
	if (potential == 0) filename += "rho";
	else if (potential == 1) filename += "m";
	std::string s = std::to_string(range); s.erase(s.size()-4,s.size());
	filename += s;
	if (sticky == 0) filename += "LOW.txt";
	else if (sticky == 1) filename += "MED.txt";
	else if (sticky == 2) filename += "HIGH.txt";

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
		//std::cout << clusters[entry] << "\n";
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

	ofile.close();
	
}

bool testSame(int N, double* clusters, int c1, int c2, double& distance) {
	//check if two clusters are the same just by comparing list of particle distances

	//set tolerance for sorted distances
	double tol = 1e-4;

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

bool testSame(int N, column_vector c1, column_vector c2, double& distance) {
	//check if two clusters are the same just by comparing list of particle distances

	//set tolerance for sorted distances
	double tol = 1e-2;

	//make particle arrays
	double* particles1 = new double[N*DIMENSION];
	double* particles2 = new double[N*DIMENSION];
	c2p(c1, particles1, N); c2p(c2, particles2, N);

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

int countNegs(column_vector eig, int N, std::vector<int>& indices) {
	//count the number of negative eigenvalues - return their indices

	int entries = N*DIMENSION;
	int negs = 0;
	double tol = 1e-8;

	for (int i = 0; i < entries; i++) {
		double eVal = eig(i);
		if (eVal < 0 && abs(eVal) > tol) {
			negs++;
			indices.push_back(i);
			//std::cout << eVal << "\n";
		}
	}

	return negs;
}

bool isMin(eigenvalue_decomposition<matrix<double>> Hd, int N, std::vector<int>& indices) {
	//check if the hessian is positive definite

	//get the eigenvalues
	column_vector e = real(Hd.get_eigenvalues());

	//check for negative eigenvalues
	int negs = countNegs(e, N, indices);

	//return false if negative eigenvalues are present
	if (negs > 0) {
		return false;
	}

	//no negative eigenvalues, this is a min
	return true;
}

void reOpt(column_vector& X, double range, double E, int N, std::vector<int> indices,
					 Parameters p, eigenvalue_decomposition<matrix<double>> Hd) {
	//re-optimize the cluster that reached a saddle

	//get the index of eigenvector with negative eigenvalue
	int index = indices[0]; 

	//construct the matrix of eigenvectors, grab the one in index, zero out rigid motion dof
	column_vector V = real(Hd.get_v());
	column_vector v = colm(V, index);
	v(0) = 0; v(1) = 0; v(2) = 0; v(4) = 0; v(5) = 0; v(8) = 0;

	//make new trial minima - left and right displacement
	column_vector Yleft(N*DIMENSION); column_vector Yright(N*DIMENSION); 


	//do re-opt left
	double p_step = 1e-10; //step to take in direction of eigenvector
	while (p_step < 1) {
		Yleft = X - p_step * v;

		//perform minimization
		find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
           objective_delta_stop_strategy(1e-14), //.be_verbose(), // Stop when the change in rosen() is less than 1e-7
           [&p](const column_vector& a) {return p.getU(a);}, 
           [&p](const column_vector& a) {return p.getGrad(a);}, 
           Yleft, -10000);

		//compute the hessian
		matrix<double> Hleft; Hleft = zeros_matrix<double>(DIMENSION*N,DIMENSION*N);
		hessMorse(Yleft, range, E, N, Hleft);

		//eigenvalue decomposition of the hessian
		eigenvalue_decomposition<matrix<double>> Hdleft(Hleft); 

		//check for minimum
		std::vector<int> neg_index;  //indices of negative eigenvalues
		if (isMin(Hdleft, N, neg_index)) {
			//printf("Escape at step %f\n", p_step);
			break;
		}
		else {
			p_step *= 5;
		}
	}

	if (p_step > 1) {
		std::cout << "Left reopt failed\n";
	}

	//do re-opt left
	p_step = 1e-10; //step to take in direction of eigenvector
	while (p_step < 1) {
		Yright = X + p_step * v;

		//perform minimization
		find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
           objective_delta_stop_strategy(1e-14), //.be_verbose(), // Stop when the change in rosen() is less than 1e-7
           [&p](const column_vector& a) {return p.getU(a);}, 
           [&p](const column_vector& a) {return p.getGrad(a);}, 
           Yright, -10000);

		//compute the hessian
		matrix<double> Hright; Hright = zeros_matrix<double>(DIMENSION*N,DIMENSION*N);
		hessMorse(Yright, range, E, N, Hright);

		//eigenvalue decomposition of the hessian
		eigenvalue_decomposition<matrix<double>> Hdright(Hright); 

		//check for minimum
		std::vector<int> neg_index;  //indices of negative eigenvalues
		if (isMin(Hdright, N, neg_index)) {
			//printf("Escape at step %f\n", p_step);
			break;
		}
		else {
			p_step *= 5;
		}
	}

	if (p_step > 1) {
		std::cout << "Right reopt failed\n";
	}

	double d;
	if (testSame(N, Yleft, Yright, d)) {
		X = Yleft;

		//std::cout << X << "\n";
	}
	else {
		std::cout << d << "\n";
		//abort();
	}




}

void descent(int N, int num_clusters, double* clusters, int sticky, int potential) {
	//perform the descent with given parameters 

	double range = 50;               //set the initial range 
	double E0 = getEnergy0(sticky);  //energy at range 50
	double E = E0;                   //energy at arbitrary range
	double kappa, kappaD;            //kappa and its derivative
	stickyF(E0, range, 1, 0, kappa, kappaD); //get initial kappa
	double end = 1.0;

	range = range - STEP;
	while (range > end) {
		//get the energy from sticky parameter value 
		E = stickyNewton(E0, range, kappa, 1); 

		//construct the parameters object 
		Parameters p = Parameters(N, potential, range, E);

		//loop over the clusters 
		for (int c = 0; c < 1; c++) {
			//get the previous cluster
			column_vector X(DIMENSION*N); 
			getCluster(N, clusters, c, X);

			//perform minimization
			find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
             objective_delta_stop_strategy(1e-14), //.be_verbose(), // Stop when the change in rosen() is less than 1e-7
             [&p](const column_vector& a) {return p.getU(a);}, 
             [&p](const column_vector& a) {return p.getGrad(a);}, 
             X, -10000);

			//compute the hessian
			matrix<double> H; H = zeros_matrix<double>(DIMENSION*N,DIMENSION*N);
			hessMorse(X, range, E, N, H);

			//eigenvalue decomposition of the hessian
			eigenvalue_decomposition<matrix<double>> Hd(H); 

			//check for minimum
			std::vector<int> neg_index;  //indices of negative eigenvalues
			if (!isMin(Hd, N, neg_index)) {
				//not a minimum, do a re-optimization
				printf("Cluster %d needs to be re-optimized at range %f\n", c, range);
				reOpt(X, range, E, N, neg_index, p, Hd);
			}

			//store the new minimum
			storeCluster(N, X, c, clusters);
		}

		//output the cluters to a file
		printClusters(N, num_clusters, clusters, range, sticky, potential);

		//test if same 
		//double d;
		//std::cout << range << ' ' << testSame(N, clusters, 0, 1, d) << ' ' << d << "\n";

		//tell user the progress
		if (abs(range - (int)range) < 1e-6)  {
			printf("Clusters at %f have been stored\n", range);
		}

		//reduce the range parameter
		range -= STEP;
		
	}
}

void parallelDescent(int N, int num_clusters, double* clusters, int sticky, int potential) {
	//perform the descent with given parameters - omp parallel

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
		#pragma omp parallel for 
		for (int c = 0; c < num_clusters; c++) {
			//get the previous cluster
			column_vector X(DIMENSION*N); 
			getCluster(N, clusters, c, X);

			//perform minimization
			find_min(bfgs_search_strategy(),  // Use BFGS search algorithm
             objective_delta_stop_strategy(1e-14), //.be_verbose(), // Stop when the change in rosen() is less than 1e-7
             [&p](const column_vector& a) {return p.getU(a);}, 
             [&p](const column_vector& a) {return p.getGrad(a);}, 
             X, -10000);

			//store the new minimum
			storeCluster(N, X, c, clusters);
		}

		//output the cluters to a file
		printClusters(N, num_clusters, clusters, range, sticky, potential);

		//test if same 
		//double d;
		//std::cout << range << ' ' << testSame(N, clusters, 0, 1, d) << ' ' << d << "\n";

		//tell user the progress
		if (abs(range - (int)range) < 1e-6)  {
			printf("Clusters at %f have been stored\n", range);
		}

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
    std::cout  << X << "\n";
    std::cout  << Y << "\n";
}

/****************************************************************************/

/* Functions for finding merges */

/****************************************************************************/


bool isMerged(int cluster, std::vector<int> merged) {
	//determine if cluster has already merged
	std::vector<int>::iterator it = std::find(merged.begin(), merged.end(), cluster);
	if (it != merged.end()) {
		return true;
	}
	return false;
}

void rMinFill(int N, int num_clusters, double* clusters, double* rMins, double rho) {
	//fill in rMin for all current clusters

	for (int i = 0; i < num_clusters; i++) {
		rMins[i] = rMin(N, i, clusters, rho);
	}

}

void findMerges(int N, int num_clusters, int sticky, int potential) {
	//determine when clusters merge

	//parameters
	double range = 50;
	double* clusters = new double[DIMENSION*N*num_clusters];
	std::vector<int> merged; //merged cluster storage
	double end = 1;
	double* rMins = new double[num_clusters]; 
	double distance;

	//output file - list of merges and ranges
	std::string filename = "../graphviz/n" + std::to_string(N);
	if (potential == 0) filename += "rho";
	else if (potential == 1) filename += "m";
	if (sticky == 0) filename += "LOW";
	else if (sticky == 1) filename += "MED";
	else if (sticky == 2) filename += "HIGH";
	filename += "merges.txt";

	std::ofstream ofile;
	ofile.open(filename);


	range = range - STEP;
	while (range > end) {
		getFinite(N, sticky, potential, range, clusters);

		for (int i = 0; i < num_clusters; i++) {
			if (!isMerged(i, merged)) {
				for (int j = i+1; j < num_clusters; j++) {
					if (!isMerged(j, merged)) {
						bool same = testSame(N, clusters, i, j, distance);
						if (same) {
							//check rMins for smooth cluster
							double r1 = rMins[i]; double r2 = rMins[j];
							//std::cout << r1 << ' ' << r2 << "\n";
								printf("Cluster %d (rmin = %f) and cluster %d (rmin = %f) merge at range %f\n", i,r1,j,r2,range);
								ofile << i << ' ' << j << ' ' << std::fixed << std::setprecision(4) << range << "\n";
								merged.push_back(j);
						}
					}
				}
			}

		}
		rMinFill(N, num_clusters, clusters, rMins, range);
		range -= STEP;
	}



	//free memory
	delete []clusters; delete []rMins;

	ofile.close();

}

double rMin(int N, int cNum, double* clusters, double rho) {
	//return the rMin value for a given cluster

	//set the tolerance for existing bonds as rough fn of rho
	double tol;
	if (rho > 30) tol = 1.01;
	if (rho <= 30 && rho > 3.5) tol = 1.03;
	if (rho <= 3.5) tol = 1.42;

	//store clusters in column vectors
	column_vector cluster(DIMENSION*N);
	getCluster(N, clusters, cNum, cluster); 

	//make particle arrays
	double* particles = new double[N*DIMENSION];
	c2p(cluster, particles, N); 

	//make vectors with distances
	std::vector<double> dist;
	double* Z = new double[DIMENSION];
	for (int i = 0; i < N; i++) {
		for (int j = i+1; j < N; j++) {
			dist.push_back(euDist(particles, i, j, N, Z)); 
		}
	}

	//find the minimum entry of the vector, greater than 1
	double rmin = 1000;
	for (int i = 0; i < dist.size(); i++) {
		if (dist[i] < rmin && dist[i] > tol) {
			rmin = dist[i];
		}
	}

	if (rmin == 1000) rmin = 1;


	//free memory
	delete []particles; delete []Z;

	//return rmin
	return rmin;
}




