#pragma once
#include <vector>
#include "string.h"
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>



void printCluster(std::ofstream& out_str, int index, std::string* range, int potential);
void printClusterID(std::ofstream& out_str, int index, std::string* range, int potential);

void sameRank(std::ofstream& out_str, std::vector<std::string> states);

void makeEdge(std::ofstream& out_str, int source, int target);
void makeInvisible(std::ofstream& out_str, std::string source, std::string target);

void makeGraph(int N, int num_clusters, int sticky, int potential);

void printCluster(std::ofstream& out_str, int index, std::string* range, std::string r2);

void makeGraph(int N, int num_clusters, int sticky); 