#pragma once
#include <memory>
#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <chrono>
#include <math.h>
#include <algorithm>
#include <string>
#include <iostream>

#define NIL -1
using namespace std;
using Item = int;
// using Transaction = set<Item>;
// using TransformedPrefixPath = std::pair<std::vector<Item>, uint64_t>;
// using TransformedPrefixPath = pair<Transaction, uint64_t>;
// using Pattern = std::pair<std::set<Item>, uint64_t>;
// using TransSet = vector<int>;
// using TransformedPrefixPath = std::pair<std::vector<Item>, TransSet>;
// using Pattern = pair<set<Item>, TransSet>;
using Edge = pair<int, int>;

using EdgesSet = vector<Edge>;
// using Clique = vector<int>;
using Frac = pair<int, int>;
// #ifndef MIN
//   #define MIN(a,b) ((a)<(b)?(a):(b))
// #endif

// #ifndef MAX
//   #define MAX(a,b) ((a)>(b)?(a):(b))
// #endif
int gcd(int num, int deno);
// 自定义比较函数，使 map 按键降序排列


struct Weight 
{
    int numerator;
    int denominator;
    Weight()
    {

    }
    Weight(int num, int deno)
    {
        int c = gcd(num, deno);
        numerator = num / c;
        denominator = deno / c;
        if (denominator < 0) {
            numerator = -numerator;
            denominator = -denominator;
        }
    }
    bool operator < (const Weight& b) const 
    {
        return numerator * b.denominator - denominator * b.numerator < 0;
    }
    bool operator < (const int& b) const
    {
        return numerator - denominator * b < 0;
    }
    bool operator <= (const Weight& b) const 
    {
        return numerator * b.denominator - denominator * b.numerator <= 0;
    }
    bool operator > (const Weight& b) const
    {
        return numerator * b.denominator - denominator * b.numerator > 0;
    }
    bool operator > (const int& b) const
    {
        return numerator - denominator * b > 0;
    }
    bool operator >= (const Weight& b) const
    {
        return numerator * b.denominator - denominator * b.numerator >= 0;
    }
    bool operator == (const Weight& b) const
    {
        return numerator == b.numerator && denominator == b.denominator;
    }
    bool operator != (const Weight& b) const
    {
        return !(numerator == b.numerator && denominator == b.denominator);
    }
    Weight operator * (const Weight& b) const
    {
        return Weight(numerator * b.numerator, denominator * b.denominator);
    }
    Weight operator * (const int& b) const
    {
        return Weight(numerator * b, denominator);
    }
    Weight operator + (const Weight& b) const
    {
        return Weight(numerator * b.denominator + b.numerator * denominator, denominator * b.denominator);
    }
    Weight operator + (const int& b) const
    {
        return Weight(numerator + b * denominator, denominator);
    }
    Weight operator - (const int& b) const
    {
        return Weight(numerator - b * denominator, denominator);
    }
    Weight operator - () const
    {
        return Weight(-numerator, denominator);
    }
    Weight operator / (const Weight& b) const
    {
        return Weight(numerator * b.denominator, denominator * b.numerator);
    }
    Weight operator / (const int& b) const
    {
        return Weight(numerator, denominator * b);
    }
    // Weight Dist() const // ��1 - a/b
    // {
    //     return Weight(denominator - numerator, denominator);
    // }
    string toString() const
    {
        return to_string(numerator) + "/" + to_string(denominator);
    }
};
const Weight WEIGHT_ONE(1, 1);
const Weight WEIGHT_ZERO(0, 1);
const Weight WEIGHT_MAX(INT16_MAX, 1);

struct CompareKeys
{
    bool operator()(const Weight &a, const Weight &b) const
    {
        return a > b;
    }
};
struct ETNode
{
    Weight score;
    set<Edge> edgeSet;
    vector<shared_ptr<ETNode>> childList;
    // unordered_map<int, set<int>> kwMap; // key:keyword; value:nodes

    ETNode(Weight score) : score(score) {}
};

struct UNode
{
    Edge value;
    shared_ptr<UNode> parent = nullptr;
    int rank = -1;
    Edge represent; // this variable is used for updating our tree index

    // Constructor
    UNode(Edge value) : value(value) {}
};

void makeSet(std::shared_ptr<UNode> &x);

std::shared_ptr<UNode> &find(std::shared_ptr<UNode> &x);

void unite(std::shared_ptr<UNode> &x, std::shared_ptr<UNode> &y);//union

void update_sup(unordered_map<int, unordered_map<int, int>> &graph, Edge &added_edge);
void MaintainKTruss_DeleteEdge(unordered_map<int, unordered_map<int, int>> &graph, vector<Edge> &temp_sup_less_edges, int k);

// unordered_map<int, pair<bool, pair<Weight, vector<int>>>> LoadAtt(string dataset, vector<int> posKws, vector<int> negKws);
void MaintainKTruss_DeleteEdge(unordered_map<int, unordered_map<int, int>> &graph, vector<Edge> &sup_less_edges, vector<Edge> &temp_sup_less_edges, int k);
Edge GetEdge(int u, int v);
struct VertexInfo
{
    int vertex;
    unordered_map<int, int> neighbors; // neighbor and support
    vector<Edge> changed_edges;
    vector<pair<Edge, int>> sup_less_edges;
};
extern  double max_time;