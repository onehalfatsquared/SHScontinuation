#include "dlib/optimization.h"
#include "continuation.h"




int main() {

	double* clusters = new double[2*DIMENSION*6];
	getSHS(6,2,clusters);
	testCV(clusters);

	delete []clusters;



	return 0;
}