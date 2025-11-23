#pragma once
#include "KTruss.h"
#include "Index.h"
#include "Utility.h"
#include <utility>
#include <stack>
struct Node {
    int vertex;   // 当前节点
    double distance; // 到该节点的距离

    // 定义优先队列的比较规则（距离小的优先）
    bool operator>(const Node& other) const {
        return distance > other.distance;
    }
};
void printf_results(DataGraph &dg, map<Weight, vector<vector<int>>, greater<Weight>> &result, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att);
void DFS(const unordered_map<int, unordered_map<int, int>> &kct, int start, unordered_set<int> &visited);
bool isCntKTruss(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct, int k);
Weight NSize(const vector<int> &M, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att);
VertexInfo MaintainTruss_DeleteVertex(int u, unordered_map<int, unordered_map<int, int>> &kct, int k, set<int> &del_vertex);
void recoverTruss(VertexInfo &vertexInfo, unordered_map<int, unordered_map<int, int>> &kct);
void max_compare(vector<int> cand, Weight &score, map<Weight, vector<vector<int>>, greater<Weight>> &result);
vector<vector<int>> FindConnectedComponents(const unordered_map<int, unordered_map<int, int>> &kct, vector<int> &cc);
bool CheckNsize(const vector<int> &M, const vector<int> &C, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, Weight sthd1, Weight sthd2);
map<Weight, vector<vector<int>>, greater<Weight>> Enumerate(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerate(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);
map<Weight, vector<vector<int>>, greater<Weight>> NoSpecial(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerateIncG(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerate1(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerate2(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);
map<Weight, vector<vector<int>>, greater<Weight>> Index_Enumerate(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);