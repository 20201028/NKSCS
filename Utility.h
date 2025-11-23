#include <sstream>
#include <fstream>
#include <iostream>
// #include "Define.h"
#include "DataGraph.h"
unordered_map<int, pair<bool, pair<Weight, vector<int>>>> LoadAtt(int method, string AttFileName, vector<int> posKws, vector<int> negKws);
DataGraph LoadGraph(string dataset, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &id2att);
Weight cal_score(Weight u_score, Weight v_score, Weight u_v_sim);
vector<vector<int>> FindConnectedComponents(const unordered_map<int, unordered_map<int, int>> &kct);
void DFS(const unordered_map<int, unordered_map<int, int>> &kct, int start, unordered_set<int> &visited, vector<int> &component);
void printGraphAsEdgeList(const std::unordered_map<int, std::unordered_map<int, int>> &graph);
unordered_map<int, vector<int>> LoadAtt(string dataset);
void calculate_score_edges(unordered_map<int, unordered_map<int, int>> &kct,
                           unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
                           map<Edge, Weight> &edge_sim,
                           map<Weight, vector<pair<int, int>>> &score_edges);