#include "graph.h"








void printCluster(std::ofstream& out_str, int index, std::string* range, int potential) {
	//print a node creation in graphviz

	std::string rangeChar;
	if (potential == 0) {
		rangeChar = "p";
		out_str << "\"" + std::to_string(index) + rangeChar + range[index] +"\" [label= \"" 
		+ "&#961;=" + range[index] +
		+ "\" , shape=circle  , width=0.1, regular=1, style=filled, fillcolor=white]; \n";
	}
	else if (potential == 1) {
		rangeChar = "m";
		out_str << "\"" + std::to_string(index) + rangeChar + range[index] +"\" [label= \"" 
		+ "m=" + range[index] +
		+ "\" , shape=circle  , width=0.1, regular=1, style=filled, fillcolor=white]; \n";
	}
	
	
} 


void sameRank(std::ofstream& out_str, std::vector<std::string> states) {
	//print a rank = same statement for all elements of states

	out_str << "{rank = same; ";
	for (int i = 0; i < states.size(); i++) {
		out_str << "\"" + states[i] + "\";";
	}
	out_str << "}\n";
}

void makeEdge(std::ofstream& out_str, std::string source, std::string target) {
	//draw an edge from source to target - no labels

	out_str << "\"" + source + "\" -> \"" + 
	target + "\" \n";
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

	//declare a file to write output to
	std::string out;
	out = "graph.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "digraph merge {\n nodesep = 1.0; ranksep = 2; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//create an array of merge values
	std::string* mergeVals = new std::string[num_clusters];
	for (int cluster = 0; cluster < num_clusters; cluster++) mergeVals[cluster] = "50.00";

	//character for range, rho or m
	std::string rangeChar;
	if (potential == 0) rangeChar = "p";
	else if (potential == 1) rangeChar = "m";

	//create vector for same rank output - create initial clusters at same rank
	std::vector<std::string> ranks; 
	for (int cluster = 0; cluster < num_clusters; cluster++) {
		printCluster(out_str, cluster, mergeVals, potential);
		ranks.push_back(std::to_string(cluster)+rangeChar+mergeVals[cluster]);
	}
	sameRank(out_str, ranks);

	//clear out ranks and set check for setting same ranks
	ranks.clear();
	std::string prev_range = "50.00";

	//go through the file 
	int i,j; double rangeD;
	while (in_str >> i) {
		//get the data in the current line
		in_str >> j; in_str >> rangeD;
		std::string range = std::to_string(rangeD); range.erase(range.size()-4,range.size());

		//check if this node already exists
		std::string current_range = range;
		std::string old_range = mergeVals[i]; 
		if (old_range != current_range) { //need to make new node for this range
			mergeVals[i] = current_range;
			printCluster(out_str, i, mergeVals, potential);
			makeEdge(out_str, std::to_string(i)+rangeChar+old_range, 
							 std::to_string(i)+rangeChar+current_range);
		}
		//if here, min index node exists, connect j to it
		old_range = mergeVals[j]; 
		mergeVals[j] = current_range;
		makeEdge(out_str, std::to_string(j)+rangeChar+old_range, 
						 std::to_string(i)+rangeChar+current_range);
		ranks.push_back(std::to_string(i)+rangeChar+current_range);
		sameRank(out_str, ranks); ranks.clear();

	}

	out_str << "}";
	delete []mergeVals;
}


void printCluster(std::ofstream& out_str, int index, std::string* range, std::string r2) {
	//print a node creation in graphviz - both ranges


	out_str << "\"" + std::to_string(index) + "p" + range[index] +"\" [label= \"" 
	+ "&#961;=" + range[index] + "\\n m=" + r2 +
	+ "\" , shape=circle  , width=0.1, regular=1, style=filled, fillcolor=white]; \n";
	
	
} 


void makeGraph(int N, int num_clusters, int sticky) {
	//make the graphviz graph - if same structure use both potentials

	//get the filename
	std::string filename1 = "../graphviz/n" + std::to_string(N);
	std::string filename2 = "../graphviz/n" + std::to_string(N);
	filename1 += "rho";
	filename2 += "m";
	if (sticky == 0) {
		filename1 += "LOW"; filename2 += "LOW";
	}
	else if (sticky == 1) {
		filename1 += "MED"; filename2 += "MED";
	}
	else if (sticky == 2) {
		filename1 += "HIGH"; filename2 += "HIGH";
	}
	filename1 += "merges.txt"; filename2 += "merges.txt";

	//open the files
	std::ifstream in_strM(filename1);
	std::ifstream in_strL(filename2);

	//check if the file can be opened
	if (!in_strM) {
		fprintf(stderr, "Cannot open file %s\n", filename1.c_str());
		return;
	}
	if (!in_strL) {
		fprintf(stderr, "Cannot open file %s\n", filename2.c_str());
		return;
	}

	//declare a file to write output to
	std::string out;
	out = "graph.txt";
	std::ofstream out_str(out);

	//write the graphviz header
	out_str << "digraph merge {\n nodesep = 1.0; ranksep = 2; \n";
	out_str << "edge [ fontcolor=red, fontsize=48];\n";

	//create an array of merge values
	std::string* mergeVals = new std::string[num_clusters];
	for (int cluster = 0; cluster < num_clusters; cluster++) mergeVals[cluster] = "50.00";

	//character for range, rho or m
	std::string rangeChar;
	rangeChar = "p";

	//create vector for same rank output - create initial clusters at same rank
	std::vector<std::string> ranks; 
	for (int cluster = 0; cluster < num_clusters; cluster++) {
		printCluster(out_str, cluster, mergeVals, "50.00");
		ranks.push_back(std::to_string(cluster)+rangeChar+mergeVals[cluster]);
	}
	sameRank(out_str, ranks);

	//clear out ranks and set check for setting same ranks
	ranks.clear();
	std::string prev_range = "50.00";

	//go through the file 
	int i,j; double rangeD;
	int i2,j2; double rangeD2;
	while (in_strM >> i) {
		//get the data in the current line - Morse
		in_strM >> j; in_strM >> rangeD;
		//get the data in the current line - LJ
		in_strL >> i2; in_strL >> j2; in_strL >> rangeD2;
		std::string range = std::to_string(rangeD); range.erase(range.size()-4,range.size());
		std::string range2 = std::to_string(rangeD2); range2.erase(range2.size()-4,range2.size());

		//check if this node already exists
		std::string current_range = range;
		std::string old_range = mergeVals[i]; 
		if (old_range != current_range) { //need to make new node for this range
			mergeVals[i] = current_range;
			printCluster(out_str, i, mergeVals, range2);
			makeEdge(out_str, std::to_string(i)+rangeChar+old_range, 
							 std::to_string(i)+rangeChar+current_range);
		}
		//if here, min index node exists, connect j to it
		old_range = mergeVals[j]; 
		mergeVals[j] = current_range;
		makeEdge(out_str, std::to_string(j)+rangeChar+old_range, 
						 std::to_string(i)+rangeChar+current_range);
		ranks.push_back(std::to_string(i)+rangeChar+current_range);
		sameRank(out_str, ranks); ranks.clear();

	}

	out_str << "}";
	delete []mergeVals;
}


