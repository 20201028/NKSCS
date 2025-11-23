#pragma once

#include <vector>
#include <string>
// #include <fstream>
#include <unordered_map>
#include <map>
// #include <iostream>
#include <set>
#include <list>
#include <queue>
#include "Define.h"

using namespace std;

// struct SimG
// {
// 	int n, m;
// 	vector<int> ei, ej;
// 	vector<double> w;
// };
class DataGraph
{
public:
	// unordered_map<int, int> graph;

	int vNum,eNum;
	vector<vector<int>> AdjList;
	// vector<int> verSort;
	unordered_map<int, int> id2seq;		 // map the vertex from ordinary dataset id to AdjList sequence. id of a node to line index of this node in adjacency list
	unordered_map<int, int> seq2id;		 // map the sequence back to the id for presentation convenience
	unordered_map<int, vector<int>> id2att; // seq and original attribute

	DataGraph();

	void addEdgeNoMatinC(int src, int dst);
	
};
