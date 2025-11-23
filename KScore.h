#pragma once
#include "Define.h"
#include "Utility.h"
class KScoreTrussDecomposition
{
public:
    map<Edge, Weight> edge2score;
    vector<int> bin;
    KScoreTrussDecomposition();
    // KScoreTrussDecomposition(unordered_map<int, unordered_map<int, int>> &graph, vector<pair<Edge, Weight>> &graph_sim, int k);
    KScoreTrussDecomposition(unordered_map<int, unordered_map<int, int>> graph, map<Edge, Weight> graph_sim, int k);
    KScoreTrussDecomposition(unordered_map<int, unordered_map<int, int>> &graph, map<Weight, vector<Edge>> &sim_graph, map<Edge, Weight> &graph_sim, int k, chrono::high_resolution_clock::time_point startTime);
};