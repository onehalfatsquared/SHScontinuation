
#include "continuation.h"



int main(int argc, char* argv[]) {

	//handle input
	if (argc != 6) {
		fprintf(stderr, "Usage: <Particles> <Kappa1> <kappa2> <pot1> <pot2> %s\n", argv[0]);
		return 1;
	}
	int N = atoi(argv[1]);         //number of particles
	int sticky1 = atoi(argv[2]);   //sticky parameter (0,1,2) -> (L,M,H)
	int sticky2 = atoi(argv[3]);   //(0,1) ->(Morse, LJ)
	int potential1 = atoi(argv[4]);   //sticky parameter (0,1,2) -> (L,M,H)
	int potential2 = atoi(argv[5]);   //(0,1) ->(Morse, LJ)

	int num_clusters = getNumClusters(N);

	getEditDistance(N, sticky1, potential1, sticky2, potential2);
	//getScatterData(N);



	return 0;
}