#pragma once
#include <algorithm>
#include <queue>
#include "DataGraph.h"
// #include "Define.h"
#define INITIAL -1
// bool compareDegree(const int &a, const int &b);
class TrussDecomposition
{
public:
    // DataGraph *datagraph;
    int kMax;
    unordered_map<int, unordered_map<int, int>> kTruss;
    unordered_map<int, unordered_map<int, int>> edge_truss;
    // map<Edge, int> trussd;
    map<int, vector<Edge>, greater<int>> k2edge;
    // map<Edge, int> sup;
    // vector<int> degree;
    // vector<int> total_order;
    // vector<int> order_pointer;
    vector<int> bin;
    TrussDecomposition();
    TrussDecomposition(DataGraph datagraph);
    // TrussDecomposition(unordered_map<int, unordered_map<int, int>> &graph, int k);
    TrussDecomposition(DataGraph datagraph, int k, chrono::high_resolution_clock::time_point startTime);
    void UpdateTrussofInsert2(vector<Edge> &, Edge &e, DataGraph *relSubG);
    // void estimateTruss(Edge &e, DataGraph *relSubG);
    void MaintainSupport_AddEdge(DataGraph *datagraph, vector<set<Edge>> &nonDecSup, Edge e);
    // bool compareDegree(const int, const int);
};

// KTruss::KTruss(DataGraph &dg) : datagraph(dg) {}
// vector<int> KTruss::degree;
// descending sort as trussness
// bool KTruss::compareDegree(const int a, const int b)
// {
//     return degree[a] > degree[b] || (degree[a] == degree[b] && a > b);
// }
