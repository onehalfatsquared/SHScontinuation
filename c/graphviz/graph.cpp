#include "graph.h"






/*

void printCluster(std::ofstream& out_str, int index, double range) {
	//print a node creation in graphviz

	out_str << "\"" + std::to_string(index) + "\" [label= \"" + std::to_string(index)
	+ "\" , shape=circle, width = 2, regular = 1, style = filled, fillcolor=white]; \n";

	%dp%d" [label="%d \\n&#961;=%d"  ,  shape=circle  , width=0.1, regular=1,style=filled,fillcolor=white]
	 ;\n',j,mergeIndex(j),j,mergeIndex(j));"
	
} 
*/

void sameRank(std::ofstream& out_str, std::vector<int> states) {
	//print a rank = same statement for all elements of states

	out_str << "{rank = same; ";
	for (int i = 0; i < states.size(); i++) {
		out_str << "\"" + std::to_string(states[i]) + "\";";
	}
	out_str << "}\n";
}

void makeEdgeClean(std::ofstream& out_str, int source, int target, double edgeWidth) {
	//draw an edge from source to target - no labels

	out_str << "\"" + std::to_string(source) + "\" -- \"" + 
	std::to_string(target) + "\" [penwidth = " + std::to_string(edgeWidth) + "]\n";
}

void makeGraph(int N, int num_clusters, int sticky, int potential) {
	//make the graphviz graph

	//get the filename
	std::string filename = "../graphviz/n" + std::to_string(N);
	if (potential == 0) filename += "rho";
	else if (potential == 1) filename += "m";
	if (sticky == 0) filename += "LOW";
	else if (sticky == 1) filename += "MED";
	else if (sticky == 2) filename += "HIGH";
	filename += "merges.txt";

	//open the file
	std::ifstream in_str(filename);

	//check if the file can be opened
	if (!in_str) {
		fprintf(stderr, "Cannot open file %s\n", filename.c_str());
		return;
	}

	//create an array of merge values
	std::string* mergeVals = new std::string[num_clusters];
	for (int cluster = 0; cluster < num_clusters; cluster++) mergeVals[cluster] = "SHS";

	//create the initial shs clusters nodes



	//make the top row
	for (int top = 0; top < num_clusters; top++) {

	}

	//go through the file 
	int i,j; double range;
	while (in_str >> i) {
		in_str >> j; in_str >> range;


	}
}