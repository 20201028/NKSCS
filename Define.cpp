#include "Define.h"
// void makeSet(shared_ptr<UNode> &x)
// {
// 	x->parent = x;
// 	x->rank = 0;
// }
// shared_ptr<UNode> find(shared_ptr<UNode> &x)
// {
// 	if (x->parent != x)
// 	{
// 		x->parent = find(x->parent); // make the parent of x point to the root
// 	}
// 	return x->parent;
// }
// void unite1(shared_ptr<UNode> &x, shared_ptr<UNode> &y) // union
// {
// 	shared_ptr<UNode> xRoot = find(x);
// 	shared_ptr<UNode> yRoot = find(y);

// 	// Compare the underlying values rather than the pointers
// 	if (xRoot->value == yRoot->value)
// 	{
// 		return;
// 	}

// 	// x and y are not already in the same set. Merge them.
// 	if (xRoot->rank < yRoot->rank)
// 	{
// 		xRoot->parent = yRoot;
// 	}
// 	else if (xRoot->rank > yRoot->rank)
// 	{
// 		yRoot->parent = xRoot;
// 	}
// 	else
// 	{
// 		yRoot->parent = xRoot;
// 		xRoot->rank = xRoot->rank + 1;
// 	}
// }

int gcd(int num, int deno)
{
	num = abs(num);
    deno = abs(deno);
    int a = max(num, deno), b = min(num, deno), c = 1;
    while (b != 0)
    {
        c = a;
        a = b;
        b = c % b;
    }
    return a;
}
void makeSet(std::shared_ptr<UNode> &x)
{
    x->parent = x;
    x->rank = 0;
    x->represent = x->value; // this is initialized as itself for our tree index
}

std::shared_ptr<UNode> &find(std::shared_ptr<UNode> &x)
{
    if (x->parent != x)
    {
        x->parent = find(x->parent);//make the parent of x point to the root
    }
    return x->parent;
}

void unite(std::shared_ptr<UNode> &x, std::shared_ptr<UNode> &y)//union
{
    std::shared_ptr<UNode> &xRoot = find(x);
    std::shared_ptr<UNode> &yRoot = find(y);

    // Compare the underlying values rather than the pointers
    if (xRoot->value == yRoot->value)
    {
        return;
    }

    // x and y are not already in the same set. Merge them.
    if (xRoot->rank < yRoot->rank)
    {
        xRoot->parent = yRoot;
    }
    else if (xRoot->rank > yRoot->rank)
    {
        yRoot->parent = xRoot;
    }
    else
    {
        yRoot->parent = xRoot;
        xRoot->rank = xRoot->rank + 1;
    }
}
void update_sup(unordered_map<int, unordered_map<int, int>> &graph, Edge &added_edge)
{
    int u = added_edge.first;
    int v = added_edge.second;

    if (!graph.count(u))
    {
        graph.emplace(u, unordered_map<int, int>());
    }
    graph[u].emplace(v, -1);
    if (!graph.count(v))
    {
        graph.emplace(v, unordered_map<int, int>());
    }
    graph[v].emplace(u, -1);
    int sup = 0;
    for (auto &it3 : graph[u])
    {
        int w = it3.first;
        if (graph[v].count(w))
        {
            sup++;
            graph[v][w]++;
            graph[w][v]++;
            graph[u][w]++;
            graph[w][u]++;
        }
    }
    graph[u][v] = sup;
    graph[v][u] = sup;
}
void MaintainKTruss_DeleteEdge(unordered_map<int, unordered_map<int, int>> &graph, vector<Edge> &temp_sup_less_edges, int k)
{
    while (!temp_sup_less_edges.empty())
    {
        auto it = temp_sup_less_edges.back();
        int u = it.first;
        int v = it.second;
        // int sup = graph[u][v];
        temp_sup_less_edges.pop_back();
        if (graph.count(u) && graph[u].count(v))
        {
            for (auto it1 : graph[v])
            {
                int w = it1.first;
                if (graph[u].count(w))
                {
                    if (graph[u][w] == k - 2)
                    {
                        temp_sup_less_edges.push_back(pair<int, int>(u, w));
                    }
                    if (graph[v][w] == k - 2)
                    {
                        temp_sup_less_edges.push_back(pair<int, int>(v, w));
                    }
                    --graph[u][w];
                    --graph[w][u];
                    --graph[v][w];
                    --graph[w][v];

                    // if (--graph[u][w] < k - 2)
                    // {
                    //     temp_sup_less_edges.push_back(pair<int, int>(u, w));
                    // }
                    // if (--graph[v][w] < k - 2)
                    // {
                    //     temp_sup_less_edges.push_back(pair<int, int>(v, w));
                    // }
                }
            }
            graph[u].erase(v);

            graph[v].erase(u);
        }
    }
}
void MaintainKTruss_DeleteEdge(unordered_map<int, unordered_map<int, int>> &graph, vector<Edge> &sup_less_edges, vector<Edge> &temp_sup_less_edges, int k)
{
    while (!temp_sup_less_edges.empty())
    {
        auto it = temp_sup_less_edges.back();
        int u = it.first;
        int v = it.second;
        // int sup = graph[u][v];
        temp_sup_less_edges.pop_back();
        if (graph.count(u) && graph[u].count(v))
        {
            // pair<int, int> v_sup(v, sup);
            sup_less_edges.push_back(pair<int, int>(u, v));
            for (auto it1 : graph[v])
            {
                int w = it1.first;
                if (graph[u].count(w))
                {
                    if (graph[u][w] == k - 2)
                    {
                        temp_sup_less_edges.push_back(pair<int, int>(u, w));
                    }
                    if (graph[v][w] == k - 2)
                    {
                        temp_sup_less_edges.push_back(pair<int, int>(v, w));
                    }
                    --graph[u][w];
                    --graph[w][u];
                    --graph[v][w];
                    --graph[w][v];

                    // if (--graph[u][w] < k - 2)
                    // {
                    //     temp_sup_less_edges.push_back(pair<int, int>(u, w));
                    // }
                    // if (--graph[v][w] < k - 2)
                    // {
                    //     temp_sup_less_edges.push_back(pair<int, int>(v, w));
                    // }
                }
            }
            graph[u].erase(v);

            graph[v].erase(u);
        }
    }
}
Edge GetEdge(int u, int v)
{
    if (u < v) return make_pair(u, v);
    return make_pair(v, u);
}
double max_time = 30.0;