#include "graph.h"
#include "continuation.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 4) {
		fprintf(stderr, "Usage: <Particles> <Kappa> <potential> %s\n", argv[0]);
		return 1;
	}
	int N = atoi(argv[1]);         //number of particles
	int kapType = atoi(argv[2]);   //sticky parameter (0,1,2) -> (L,M,H)
	int potential = atoi(argv[3]);   //(0,1) ->(Morse, LJ)

	int num_clusters = getNumClusters(N);

	//makeGraph(N, num_clusters, kapType, potential);
	makeGraph(N, num_clusters, kapType);

	



	return 0;
}