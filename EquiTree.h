#pragma once
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "KTruss.h"
#include "KScore.h"
// #include "BasicIndex.h"
using namespace std;

// struct AdvancedIndex
// {
//     map<int, shared_ptr<ETNode>> kScore;

//     // shared_ptr<ETNode> root;
//     AdvancedIndex(DataGraph &);
// };
// bool compareDegree1(const int &a, const int &b, DataGraph &dg)
// {
//     return dg.AdjList[a].size() > dg.AdjList[b].size() || (dg.AdjList[a].size() == dg.AdjList[b].size() && a > b);
// }


map<int, shared_ptr<ETNode>> AdvancedIndex(DataGraph &dg);