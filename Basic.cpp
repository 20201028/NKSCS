#include "Basic.h"
#include "KScore.h"
#include <algorithm>
#include <sstream>

// unordered_map<int, vector<int>> LoadAtt(string dataset)
// {
//     unordered_map<int, vector<int>> id2att;
//     string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
//     ifstream ain(AttFileName.c_str());

//     if (!ain)
//     {
//         cout << "Fail to read " << AttFileName << "." << endl;
//         return {};
//     }

//     string aline;

//     while (getline(ain, aline))
//     {
//         if (aline.find('#') != string::npos)
//             continue;
//         string ver_s = aline.substr(0, aline.find("\t"));
//         string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
//         int ver = stoi(ver_s); // string2float
//         int start = 0, end = 0;
//         vector<int> att;
//         int w = 0;
//         while ((end = att_s.find(",", start)) != string::npos)
//         {
//             w = stoi(att_s.substr(start, end - start));
//             att.push_back(w);
//             start = end + 1;
//         }
//         w = stoi(att_s.substr(start));
//         att.push_back(w);
//         id2att.emplace(ver, att);
//     }
//     ain.close();
//     cout << "Loaded attribute successfully!" << endl;
//     return move(id2att);
// }

// DataGraph LoadGraph(string dataset, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &id2att)
// {
//     DataGraph datagraph;
//     string StructFileName = "DataGraph/" + dataset + "/" + "graph.txt";
//     ifstream sin(StructFileName.c_str());

//     if (!sin)
//     {
//         cout << "Fail to read " << StructFileName << "." << endl;
//         return datagraph;
//     }
//     string sline;
//     while (getline(sin, sline)) // 默认数据集中边不重复
//     {
//         if (sline.find('#') != string::npos)
//             continue;
//         string src_s = sline.substr(0, sline.find("\t"));
//         string dst_s = sline.substr(sline.find("\t") + 1, sline.find("\n") - sline.find("\t") - 1);
//         int src = stoi(src_s); // string2float
//         int dst = stoi(dst_s);
//         if (id2att.count(src) && id2att.count(dst))
//         {
//             datagraph.addEdgeNoMatinC(src, dst);
//         }
//     }

//     cout << "Loaded dataset successfully!" << endl;
//     return move(datagraph);
// }
// DataGraph LoadGraph(string dataset)
// {
//     DataGraph datagraph;
//     string StructFileName = "DataGraph/" + dataset + "/" + "graph.txt";
//     ifstream sin(StructFileName.c_str());

//     if (!sin)
//     {
//         cout << "Fail to read " << StructFileName << "." << endl;
//         return datagraph;
//     }
//     string sline;
//     while (getline(sin, sline)) // 默认数据集中边不重复
//     {
//         if (sline.find('#') != string::npos)
//             continue;
//         string src_s = sline.substr(0, sline.find("\t"));
//         string dst_s = sline.substr(sline.find("\t") + 1, sline.find("\n") - sline.find("\t") - 1);
//         int src = stoi(src_s); // string2float
//         int dst = stoi(dst_s);

//         datagraph.addEdgeNoMatinC(src, dst);
//     }

//     cout << "Loaded graph successfully!" << endl;

//     string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
//     ifstream ain(AttFileName.c_str());

//     if (!ain)
//     {
//         cout << "Fail to read " << AttFileName << "." << endl;
//         return datagraph;
//     }

//     string aline;

//     while (getline(ain, aline))
//     {
//         if (aline.find('#') != string::npos)
//             continue;
//         string ver_s = aline.substr(0, aline.find("\t"));
//         string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
//         int ver = stoi(ver_s); // string2float
//         if (datagraph.id2seq.count(ver))
//         {
//             int start = 0, end = 0;
//             vector<int> att;
//             int w = 0;
//             while ((end = att_s.find(",", start)) != string::npos)
//             {
//                 w = stoi(att_s.substr(start, end - start));
//                 att.push_back(w);
//                 start = end + 1;
//             }
//             w = stoi(att_s.substr(start));
//             att.push_back(w);
//             datagraph.id2att.emplace(ver, att);
//         }
//     }
//     ain.close();
//     cout << "Loaded attribute successfully!" << endl;
//     return move(datagraph);
// }
Weight NSize(const vector<int> &M, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    int n = 0;
    for (auto it : M)
    {
        if (!seq2att[it].first)
        {
            n++;
        }
    }
    return Weight(n, M.size());
}
void DFS(const unordered_map<int, unordered_map<int, int>> &kct, int start, unordered_set<int> &visited)
{
    stack<int> stack;
    stack.push(start);
    while (!stack.empty())
    {
        int node = stack.top();
        stack.pop();
        if (visited.find(node) == visited.end())
        {
            visited.insert(node);
            for (const auto &neighbor : kct.at(node))
            {
                stack.push(neighbor.first);
            }
        }
    }
}
// void DFS(const unordered_map<int, unordered_map<int, int>> &kct, int start, unordered_set<int> &visited, vector<int> &component)
// {
//     stack<int> stack;
//     stack.push(start);
//     while (!stack.empty())
//     {
//         int node = stack.top();
//         stack.pop();
//         if (visited.find(node) == visited.end())
//         {
//             visited.insert(node);
//             component.push_back(node);
//             for (const auto &neighbor : kct.at(node))
//             {
//                 stack.push(neighbor.first);
//             }
//         }
//     }
// }
// 找到所有连通分量的函数
// vector<vector<int>> FindConnectedComponents(const unordered_map<int, unordered_map<int, int>> &kct)
// {
//     vector<vector<int>> components;
//     unordered_set<int> visited;

//     for (const auto &node : kct)
//     {
//         if (visited.find(node.first) == visited.end())
//         {
//             vector<int> component;
//             DFS(kct, node.first, visited, component);
//             components.push_back(component);
//         }
//     }

//     return move(components);
// }
vector<vector<int>> FindConnectedComponents(const unordered_map<int, unordered_map<int, int>> &kct, vector<int> &cc)
{
    vector<vector<int>> components;
    unordered_set<int> visited;

    for (int v : cc)
    {
        if (visited.find(v) == visited.end())
        {
            vector<int> component;
            DFS(kct, v, visited, component);
            components.push_back(component);
        }
    }

    return move(components);
}

bool isCntKTruss(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct, int k)
{
    int scale = M.size();
    if (scale < 2)
    {
        return false; // 空图和孤立点不认为是连通的
    }
    // 检查边权重
    vector<pair<Edge, int>> sup_less_edges;
    for (int v : M)
    {
        for (const auto &neighbor : kct[v])
        {
            if (v < neighbor.first && neighbor.second < k - 2)
            {
                sup_less_edges.push_back(pair<Edge, int>(Edge(v, neighbor.first), 0));
            }
        }
    }
    // delete the edges whose sup < k-2
    for (int i = 0; i < sup_less_edges.size(); i++)
    {
        Edge e = sup_less_edges[i].first;
        int u = e.first;
        int v = e.second;
        sup_less_edges[i].second = kct[u][v];
        kct[u].erase(v);
        kct[v].erase(u);
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                if (kct[u][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(u, w), 0));
                }
                --kct[u][w];
                --kct[w][u];
                if (kct[v][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(v, w), 0));
                }
                --kct[v][w];
                --kct[w][v];
            }
        }
    }

    // 检查连通性
    unordered_set<int> visited;
    // int startNode = kct.begin()->first;

    int startNode;
    // if (M.empty())
    // {
    //     startNode = C[0];
    // }
    // else
    // {
    //     startNode = M[0];
    // }
    startNode = M[0];
    DFS(kct, startNode, visited);
    // recover
    while (!sup_less_edges.empty())
    {
        auto &it = sup_less_edges.back();

        Edge e = it.first;
        int u = e.first;
        int v = e.second;
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                ++kct[u][w];
                ++kct[w][u];
                ++kct[v][w];
                ++kct[w][v];
            }
        }
        kct[u].emplace(v, it.second);
        kct[v].emplace(u, it.second);
        sup_less_edges.pop_back();
    }

    // 检查是否所有节点都被访问
    // if (visited.size() != kct.size())
    if (visited.size() != scale)
    {
        return false; // 图不连通
    }

    return true; // 图连通且没有权重小于 k 的边
}
bool isCntKTruss2(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct, int k)
{
    int scale = M.size();
    if (scale < 2)
    {
        return false; // 空图和孤立点不认为是连通的
    }
    // 检查边权重
    vector<pair<Edge, int>> sup_less_edges;
    for (int v : M)
    {
        for (const auto &neighbor : kct[v])
        {
            if (v < neighbor.first && neighbor.second < k - 2)
            {
                sup_less_edges.push_back(pair<Edge, int>(Edge(v, neighbor.first), 0));
            }
        }
    }
    // delete the edges whose sup < k-2
    for (int i = 0; i < sup_less_edges.size(); i++)
    {
        Edge e = sup_less_edges[i].first;
        int u = e.first;
        int v = e.second;
        sup_less_edges[i].second = kct[u][v];
        kct[u].erase(v);
        kct[v].erase(u);
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                if (kct[u][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(u, w), 0));
                }
                --kct[u][w];
                --kct[w][u];
                if (kct[v][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(v, w), 0));
                }
                --kct[v][w];
                --kct[w][v];
            }
        }
    }
    bool del_edge = false;
    if (!sup_less_edges.empty())
    {
        del_edge = true;
    }
    unordered_set<int> visited;
    if (del_edge)
    { // 检查连通性

        // int startNode = kct.begin()->first;

        int startNode = M[0];
        DFS(kct, startNode, visited);
    }

    // recover
    while (!sup_less_edges.empty())
    {
        auto &it = sup_less_edges.back();

        Edge e = it.first;
        int u = e.first;
        int v = e.second;
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                ++kct[u][w];
                ++kct[w][u];
                ++kct[v][w];
                ++kct[w][v];
            }
        }
        kct[u].emplace(v, it.second);
        kct[v].emplace(u, it.second);
        sup_less_edges.pop_back();
    }

    // 检查是否所有节点都被访问
    // if (visited.size() != kct.size())
    if (del_edge)
    {
        if (visited.size() != scale)
        {
            return false; // 图不连通
        }
    }

    return true; // 图连通且没有权重小于 k 的边
}
Weight KScore(vector<int> &M, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, Weight opt)
{
    Weight min_score = WEIGHT_MAX;
    for (int i = 0; i < M.size(); i++)
    {
        int u = M[i];
        Weight u_score = seq2att[u].second.first;
        vector<int> u_att = seq2att[u].second.second;
        for (int j = i + 1; j < M.size(); j++)
        {
            int v = M[j];
            Weight v_score = seq2att[v].second.first;
            vector<int> v_att = seq2att[v].second.second;
            vector<int> com_att;
            set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
            Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
            Weight u_v = cal_score(u_score, v_score, u_v_sim);
            if (u_v < opt)
                return WEIGHT_ZERO;
            min_score = min(min_score, u_v);
        }
    }
    return min_score;
}
void max_compare(vector<int> cand, Weight &score, map<Weight, vector<vector<int>>, greater<Weight>> &result_map)
{
    if (result_map.size() == 0)
    {
        result_map.emplace(score, vector<vector<int>>(1, cand));
        // result.push_back(make_pair(score, cand));
    }
    else
    {
        auto &result = result_map.begin()->second;
        int id = -1;
        for (int i = 0; i < result.size(); i++)
        {
            // vector<int> &members = result[i].second;
            vector<int> &members = result[i];
            vector<int> com;
            set_intersection(cand.begin(), cand.end(), members.begin(), members.end(), back_inserter(com));
            if (com.size() == members.size())
            {                   //|com|<=|result[i]|
                members = cand; // cand>>result[i]
                id = i;
                // cout << "cand score: " << score << endl;
                break;
            }
            else
            {
                if (com.size() == cand.size()) // result[i]>>cand
                    break;
                else
                {
                    if (i == result.size() - 1)
                    {
                        // result.push_back(make_pair(score, cand));
                        result_map[score].push_back(cand);
                        // cout << "cand score: " << score << endl;
                    }
                }
            }
        }
        if (id != -1 && id != result.size() - 1)
        {
            auto it = result[id + 1];
            for (int i = id + 1; i < result.size();)
            {
                vector<int> &members = result[i];
                vector<int> com;
                set_intersection(cand.begin(), cand.end(), members.begin(), members.end(), back_inserter(com));
                if (com.size() == members.size())
                { //|com|<=|result[i]|
                    result.erase(move(result.begin() + i));
                    // cout << "cand score: " << score << endl;
                }
                else
                {
                    i++;
                }
            }
        }
    }
}
void record_results(vector<int> cand, Weight &score, Weight &opt, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{

    sort(cand.begin(), cand.end());
    if (score == opt)
    {
        max_compare(cand, score, result);
    }
    else
    {
        result.clear();
        // result.push_back(make_pair(score, cand));
        result.emplace(score, vector<vector<int>>(1, cand));
        // cout << "cand score: " << score << endl;
        opt = score;
    }
}
void NaiveEnum(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
               unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2 && isCntKTruss(M, C, kct, k))
        {
            Weight score = KScore(M, seq2att, opt);

            if (score > WEIGHT_ZERO)
            {
                record_results(M, score, opt, result);
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
            }
        }
        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    NaiveEnum(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
    M.pop_back(); // 恢复 M

    unordered_map<int, int> u_neighbors = kct[u];
    // local update support
    vector<Edge> changed_edges;
    set<int> u_neighbors_set;
    for (auto &v : u_neighbors)
    {
        u_neighbors_set.insert(v.first);
    }
    for (auto &v : u_neighbors)
    {
        set<int> v_neighbors_set;
        for (auto &it : kct[v.first])
        {
            v_neighbors_set.insert(it.first);
        }
        set<int> com_neighbors;
        set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
        for (int w : com_neighbors)
        {
            if (find(changed_edges.begin(), changed_edges.end(), Edge(v.first, w)) == changed_edges.end())
            {
                changed_edges.push_back(Edge(v.first, w));
                changed_edges.push_back(Edge(w, v.first));
                kct[w][v.first]--;
                kct[v.first][w]--;
            }
        }
    }
    // delete u from kct
    kct.erase(u);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].erase(u);
        }
    }

    // Shrink: 调用 NaiveEnum(M, C \ u)
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    NaiveEnum(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);

    // 恢复 C
    C.push_back(u);
    kct.emplace(u, u_neighbors);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].emplace(u, it.second);
        }
    }
    // recover support
    for (auto &e : changed_edges)
    {
        kct[e.first][e.second]++;
    }
}
pair<pair<int, unordered_map<int, int>>, unordered_map<int, vector<int>>> MaintainSupport_DeleteVertex(int u, unordered_map<int, unordered_map<int, int>> &kct)
{
    unordered_map<int, int> u_neighbors = kct[u];
    // local update support
    pair<int, unordered_map<int, int>> u_neighbors_pair(u, u_neighbors);
    unordered_map<int, vector<int>> changed_edges;
    set<int> u_neighbors_set;
    for (auto &v : u_neighbors)
    {
        u_neighbors_set.insert(v.first);
    }
    for (auto &v : u_neighbors)
    {
        set<int> v_neighbors_set;
        for (auto &it : kct[v.first])
        {
            v_neighbors_set.insert(it.first);
        }
        set<int> com_neighbors;
        set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
        if (!com_neighbors.empty())
        {
            changed_edges.emplace(v.first, vector<int>());
        }
        for (int w : com_neighbors)
        {
            changed_edges[v.first].push_back(w);
            kct[w][v.first]--;
            kct[v.first][w]--;
        }
        kct[v.first].erase(u);
        u_neighbors_set.erase(v.first);
    }
    // delete u from kct
    kct.erase(u);
    pair<pair<int, unordered_map<int, int>>, unordered_map<int, vector<int>>> ret(u_neighbors_pair, changed_edges);
    return move(ret);
}
void RecoverEdges(pair<pair<int, unordered_map<int, int>>, unordered_map<int, vector<int>>> &ret, unordered_map<int, unordered_map<int, int>> &kct)
{
    int v = ret.first.first;
    kct[v] = ret.first.second;
    for (auto it : ret.first.second)
    {
        kct[it.first][v] = it.second;
    }
    for (auto &e : ret.second)
    {
        int v = e.first;
        for (auto &w : e.second)
        {
            kct[v][w]++;
            kct[w][v]++;
        }
    }
}
void RecoverEdges(int v, unordered_map<int, int> &neighbors, vector<Edge> changed_edges, unordered_map<int, unordered_map<int, int>> &kct)
{
    kct[v] = neighbors;
    for (auto it : neighbors)
    {
        kct[it.first][v] = it.second;
    }
    for (auto &e : changed_edges)
    {
        int v = e.first;
        int w = e.second;
        kct[v][w]++;
        kct[w][v]++;
    }
}

vector<int> Shrink_C_Att(int u, Weight &opt, vector<int> &C, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    vector<int> removed_vertices;
    Weight u_score = seq2att[u].second.first;
    vector<int> u_att = seq2att[u].second.second;
    for (auto it = C.begin(); it != C.end();)
    {
        int v = *it;
        Weight v_score = seq2att[v].second.first;
        vector<int> v_att = seq2att[v].second.second;
        vector<int> com_att;
        set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
        Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
        Weight u_v = cal_score(u_score, v_score, u_v_sim);
        if (u_v < opt)
        {
            removed_vertices.push_back(v);
            it = C.erase(it); // 删除元素并更新迭代器
        }
        else
        {
            ++it; // 继续遍历下一个元素
        }
    }
    return move(removed_vertices);
}
vector<int> Shrink_C_Att(int u, Weight &opt, vector<int> &C, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
                         unordered_map<int, unordered_map<int, Weight>> &graph_sim)
{
    vector<int> removed_vertices;
    for (auto it = C.begin(); it != C.end();)
    {
        int v = *it;
        if (graph_sim[u][v] < opt)
        {
            removed_vertices.push_back(v);
            it = C.erase(it); // 删除元素并更新迭代器
        }
        else
        {
            ++it; // 继续遍历下一个元素
        }
    }
    return move(removed_vertices);
}
void NaiveEnum_Att(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                   unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k,
                   Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // result stores the result whose score is larger than opt
    Weight M_min_score = WEIGHT_MAX;
    // unordered_map<int, unordered_map<int,int>> min_score_map;
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2 && isCntKTruss(M, C, kct, k))
        {
            Weight score = KScore(M, seq2att, opt); // can further improve because the score of vertex pair is already calculated

            if (score > WEIGHT_ZERO)
            {
                record_results(M, score, opt, result);
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
            }
        }
        return;
    }
    // C\u
    int u = C.back();
    C.pop_back();

    // MUu
    M.push_back(u);

    // kscore shrink C
    // recode the min score of M, if the score is smaller than opt, then return
    vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att);
    vector<pair<pair<int, unordered_map<int, int>>, unordered_map<int, vector<int>>>> changed_edges;
    for (int v : removed_vertices)
    {
        changed_edges.push_back(MaintainSupport_DeleteVertex(v, kct));
    }

    // if(!removed_vertices.empty())
    // cout << "removed_vertices: " << removed_vertices.size() << endl;

    // ktruss

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    NaiveEnum_Att(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
    // M, C\u
    M.pop_back();
    // recover the vertices and vedges dissimilar to u
    for (int i = 0; i < removed_vertices.size(); i++)
    {
        C.push_back(removed_vertices[i]);
    }
    while (!changed_edges.empty())
    {
        pair<pair<int, unordered_map<int, int>>, unordered_map<int, vector<int>>> ret = changed_edges.back();
        RecoverEdges(ret, kct);
        changed_edges.pop_back();
    }

    // delete u
    unordered_map<int, int> u_neighbors = kct[u];
    pair<pair<int, unordered_map<int, int>>, unordered_map<int, vector<int>>> ret = MaintainSupport_DeleteVertex(u, kct);
    // Shrink: 调用 NaiveEnum(M, C \ u)
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    NaiveEnum_Att(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);

    // 恢复原状M, C
    C.push_back(u);
    RecoverEdges(ret, kct);
}

void NaiveEnum_Truss(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                     unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            // check connectivity
            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() != M.size())
            {
                return; // 图不连通
            }
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
                record_results(M, score, opt, result);
            }
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    NaiveEnum_Truss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
    M.pop_back(); // 恢复 M
    // MaintainTruss_DeleteVertex();

    // delete vertex
    unordered_map<int, int> u_neighbors = kct[u];
    // local update support
    vector<Edge> changed_edges;
    vector<pair<Edge, int>> sup_less_edges;
    set<int> u_neighbors_set;
    for (auto &v : u_neighbors)
    {
        u_neighbors_set.insert(v.first);
    }
    for (auto &v : u_neighbors)
    {
        set<int> v_neighbors_set;
        for (auto &it : kct[v.first])
        {
            v_neighbors_set.insert(it.first);
        }
        set<int> com_neighbors;
        set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
        for (int w : com_neighbors)
        {
            if (find(changed_edges.begin(), changed_edges.end(), Edge(v.first, w)) == changed_edges.end())
            {
                changed_edges.push_back(Edge(v.first, w));
                changed_edges.push_back(Edge(w, v.first));
                if (--kct[w][v.first] < (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(w, v.first), 0));
                }
                kct[v.first][w]--;
            }
        }
    }
    // delete u from kct
    kct.erase(u);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].erase(u);
        }
    }

    // delete the edges whose sup < k-2
    set<int> del_vertex;
    for (int i = 0; i < sup_less_edges.size(); i++)
    {
        Edge e = sup_less_edges[i].first;
        int u = e.first;
        int v = e.second;
        sup_less_edges[i].second = kct[u][v];
        kct[u].erase(v);
        kct[v].erase(u);
        // delete isolated vertices
        if (kct[u].empty())
            del_vertex.insert(u);
        if (kct[v].empty())
            del_vertex.insert(v);
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                if (kct[u][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(u, w), 0));
                }
                --kct[u][w];
                --kct[w][u];
                if (kct[v][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(v, w), 0));
                }
                --kct[v][w];
                --kct[w][v];
            }
        }
    }
    set<int> M_temp(M.begin(), M.end());
    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        if (M.size() + C.size() >= k)
        {
            NaiveEnum_Truss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        }
        for (int v : del_vertex)
            C.push_back(v);
    }

    // 恢复 C
    C.push_back(u);
    while (!sup_less_edges.empty())
    {
        auto &it = sup_less_edges.back();

        Edge e = it.first;
        int u = e.first;
        int v = e.second;
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                ++kct[u][w];
                ++kct[w][u];
                ++kct[v][w];
                ++kct[w][v];
            }
        }
        kct[u].emplace(v, it.second);
        kct[v].emplace(u, it.second);
        sup_less_edges.pop_back();
    }
    kct.emplace(u, u_neighbors);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].emplace(u, it.second);
        }
    }

    // recover support
    for (auto &e : changed_edges)
    {
        kct[e.first][e.second]++;
    }
}

void NaiveEnum_Cnt(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                   unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2 && isCntKTruss2(M, C, kct, k))
        {
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
                record_results(M, score, opt, result);
            }
        }
        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    NaiveEnum_Cnt(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
    M.pop_back(); // 恢复 M

    unordered_map<int, int> u_neighbors = kct[u];
    // local update support
    vector<Edge> changed_edges;
    set<int> u_neighbors_set;
    for (auto &v : u_neighbors)
    {
        u_neighbors_set.insert(v.first);
    }
    for (auto &v : u_neighbors)
    {
        set<int> v_neighbors_set;
        for (auto &it : kct[v.first])
        {
            v_neighbors_set.insert(it.first);
        }
        set<int> com_neighbors;
        set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
        for (int w : com_neighbors)
        {
            if (find(changed_edges.begin(), changed_edges.end(), Edge(v.first, w)) == changed_edges.end())
            {
                changed_edges.push_back(Edge(v.first, w));
                changed_edges.push_back(Edge(w, v.first));
                kct[w][v.first]--;
                kct[v.first][w]--;
            }
        }
    }
    // delete u from kct
    kct.erase(u);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].erase(u);
        }
    }
    // find connected compontents
    if (M.empty())
    {
        vector<vector<int>> components = FindConnectedComponents(kct, C);
        if (components.size() == 1)
        {
            NaiveEnum_Cnt(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        }
        else
        {
            for (auto &c : components)
            {
                unordered_map<int, unordered_map<int, int>> local_kct;
                for (auto v : c)
                {
                    local_kct.emplace(v, kct.at(v));
                }
                NaiveEnum_Cnt(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result);
            }
        }
    }
    else
    {

        unordered_set<int> visited;
        int startNode = M[0];
        DFS(kct, startNode, visited);
        if (visited.size() == M.size() + C.size())
        {
            NaiveEnum_Cnt(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        }
        else
        {
            // vector<int> C_temp;
            // int M_count = 0;
            // for (int v : visited)
            // {
            //     if (count(M.begin(), M.end(), v))
            //         M_count++;
            //     else
            //         C_temp.push_back(v);

            // }
            // if (M_count == M.size())
            // {
            //     unordered_map<int, unordered_map<int, int>> local_kct;
            //     for (auto v : visited)
            //     {
            //         local_kct.emplace(v, kct.at(v));
            //     }
            //     NaiveEnum_Cnt(M, C_temp, local_kct, seq2att, k, opt, sthd, result);
            // }
            bool M_count = true;
            for (auto it : M)
            {
                if (!visited.count(it))
                {
                    M_count = false;
                }
            }
            // for (int v : visited)
            // {
            //     if (count(M.begin(), M.end(), v))
            //         M_count++;
            // }
            if (M_count)
            {
                unordered_map<int, unordered_map<int, int>> delete_set;
                auto it = C.begin();
                while (it != C.end())
                {
                    if (!visited.count(*it))
                    {

                        delete_set.emplace(*it, kct.at(*it));
                        kct.erase(*it);
                        it = C.erase(it);
                    }
                    else
                        ++it;
                }
                NaiveEnum_Cnt(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
                for (auto v : delete_set)
                {
                    C.push_back(v.first);
                    kct.emplace(v.first, v.second);
                }
            }
        }
    }

    // 恢复 C
    C.push_back(u);
    kct.emplace(u, u_neighbors);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].emplace(u, it.second);
        }
    }
    // recover support
    for (auto &e : changed_edges)
    {
        kct[e.first][e.second]++;
    }
}
void NaiveEnum_CntFirst(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                        unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd, vector<vector<int>> &result)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        if (NSize(M, seq2att) >= sthd && isCntKTruss2(M, C, kct, k))
        {
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                vector<int> cand;
                for (const auto &it : M)
                {
                    cand.push_back(it);
                }
                sort(cand.begin(), cand.end());
                if (score == opt)
                {
                    for (int i = 0; i < result.size(); i++)
                    {
                        vector<int> com;
                        set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                        if (com.size() == result[i].size())
                            result[i] = cand;
                        else if (com.size() == cand.size())
                        {
                        }
                        else
                            result.push_back(cand);
                    }
                }
                else
                {
                    result.clear();
                    result.push_back(cand);
                    opt = score;
                }
            }
        }
        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);
    // find connected compontents
    if (!M.empty())
    {

        unordered_set<int> visited;
        int startNode = M[0];
        DFS(kct, startNode, visited);
        if (visited.size() == M.size() + C.size())
        {
            NaiveEnum_CntFirst(M, C, kct, seq2att, k, opt, sthd, result);
        }
        else
        {
            // vector<int> C_temp;
            // int M_count = 0;
            // for (int v : visited)
            // {
            //     if (count(M.begin(), M.end(), v))
            //         M_count++;
            //     else
            //         C_temp.push_back(v);

            // }
            // if (M_count == M.size())
            // {
            //     unordered_map<int, unordered_map<int, int>> local_kct;
            //     for (auto v : visited)
            //     {
            //         local_kct.emplace(v, kct.at(v));
            //     }
            //     NaiveEnum_Cnt(M, C_temp, local_kct, seq2att, k, opt, sthd, result);
            // }
            bool M_count = true;
            for (auto it : M)
            {
                if (!visited.count(it))
                {
                    M_count = false;
                }
            }
            // for (int v : visited)
            // {
            //     if (count(M.begin(), M.end(), v))
            //         M_count++;
            // }
            if (M_count)
            {
                unordered_map<int, unordered_map<int, int>> delete_set;
                auto it = C.begin();
                while (it != C.end())
                {
                    if (!visited.count(*it))
                    {

                        delete_set.emplace(*it, kct.at(*it));
                        kct.erase(*it);
                        it = C.erase(it);
                    }
                    else
                        ++it;
                }
                NaiveEnum_CntFirst(M, C, kct, seq2att, k, opt, sthd, result);
                for (auto v : delete_set)
                {
                    C.push_back(v.first);
                    kct.emplace(v.first, v.second);
                }
            }
        }
    }

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    M.pop_back(); // 恢复 M

    unordered_map<int, int> u_neighbors = kct[u];
    // local update support
    vector<Edge> changed_edges;
    set<int> u_neighbors_set;
    for (auto &v : u_neighbors)
    {
        u_neighbors_set.insert(v.first);
    }
    for (auto &v : u_neighbors)
    {
        set<int> v_neighbors_set;
        for (auto &it : kct[v.first])
        {
            v_neighbors_set.insert(it.first);
        }
        set<int> com_neighbors;
        set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
        for (int w : com_neighbors)
        {
            if (find(changed_edges.begin(), changed_edges.end(), Edge(v.first, w)) == changed_edges.end())
            {
                changed_edges.push_back(Edge(v.first, w));
                changed_edges.push_back(Edge(w, v.first));
                kct[w][v.first]--;
                kct[v.first][w]--;
            }
        }
    }
    // delete u from kct
    kct.erase(u);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].erase(u);
        }
    }
    NaiveEnum_CntFirst(M, C, kct, seq2att, k, opt, sthd, result);

    // 恢复 C
    C.push_back(u);
    kct.emplace(u, u_neighbors);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].emplace(u, it.second);
        }
    }
    // recover support
    for (auto &e : changed_edges)
    {
        kct[e.first][e.second]++;
    }
}
void NaiveEnum_CntTruss(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                        unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
                record_results(M, score, opt, result);
            }
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    NaiveEnum_CntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
    M.pop_back(); // 恢复 M
    // MaintainTruss_DeleteVertex();

    // delete vertex
    unordered_map<int, int> u_neighbors = kct[u];
    // local update support
    vector<Edge> changed_edges;
    vector<pair<Edge, int>> sup_less_edges;
    set<int> u_neighbors_set;
    for (auto &v : u_neighbors)
    {
        u_neighbors_set.insert(v.first);
    }
    for (auto &v : u_neighbors)
    {
        set<int> v_neighbors_set;
        for (auto &it : kct[v.first])
        {
            v_neighbors_set.insert(it.first);
        }
        set<int> com_neighbors;
        set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
        for (int w : com_neighbors)
        {
            if (find(changed_edges.begin(), changed_edges.end(), Edge(v.first, w)) == changed_edges.end())
            {
                changed_edges.push_back(Edge(v.first, w));
                changed_edges.push_back(Edge(w, v.first));
                if (--kct[w][v.first] < (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(w, v.first), 0));
                }
                kct[v.first][w]--;
            }
        }
    }
    // delete u from kct
    kct.erase(u);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].erase(u);
        }
    }

    // delete the edges whose sup < k-2
    set<int> del_vertex;
    for (int i = 0; i < sup_less_edges.size(); i++)
    {
        Edge e = sup_less_edges[i].first;
        int u = e.first;
        int v = e.second;
        sup_less_edges[i].second = kct[u][v];
        kct[u].erase(v);
        kct[v].erase(u);
        // delete isolated vertices
        if (kct[u].empty())
            del_vertex.insert(u);
        if (kct[v].empty())
            del_vertex.insert(v);
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                if (kct[u][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(u, w), 0));
                }
                --kct[u][w];
                --kct[w][u];
                if (kct[v][w] == (k - 2))
                {
                    sup_less_edges.push_back(make_pair(Edge(v, w), 0));
                }
                --kct[v][w];
                --kct[w][v];
            }
        }
    }
    set<int> M_temp(M.begin(), M.end());
    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            kct.erase(v);
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        /////////////////////////////////// find connected compontents
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_CntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
            }
            else
            {
                for (auto &c : components)
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_CntTruss(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result);
                }
            }
        }
        else
        {

            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                NaiveEnum_CntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_CntTruss(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result);
                }
            }
        }
        //////////////////////////////////////
        for (int v : del_vertex)
        {
            kct.emplace(v, unordered_map<int, int>());
            C.push_back(v);
        }
    }

    // 恢复 C
    C.push_back(u);
    while (!sup_less_edges.empty())
    {
        auto &it = sup_less_edges.back();

        Edge e = it.first;
        int u = e.first;
        int v = e.second;
        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct[v].count(w))
            {
                ++kct[u][w];
                ++kct[w][u];
                ++kct[v][w];
                ++kct[w][v];
            }
        }
        kct[u].emplace(v, it.second);
        kct[v].emplace(u, it.second);
        sup_less_edges.pop_back();
    }
    kct.emplace(u, u_neighbors);
    for (const auto &it : u_neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].emplace(u, it.second);
        }
    }

    // recover support
    for (auto &e : changed_edges)
    {
        kct[e.first][e.second]++;
    }
}

VertexInfo MaintainTruss_DeleteVertex(int u, unordered_map<int, unordered_map<int, int>> &kct, int k, set<int> &del_vertex)
{
    VertexInfo vertexInfo;
    vertexInfo.vertex = u;
    // delete vertex
    vertexInfo.neighbors = kct.at(u);
    // local update support
    set<int> u_neighbors_set;
    for (auto &v : vertexInfo.neighbors)
    {
        u_neighbors_set.insert(v.first);
    }
    for (auto &v : vertexInfo.neighbors)
    {
        set<int> v_neighbors_set;
        for (auto &it : kct.at(v.first))
        {
            v_neighbors_set.insert(it.first);
        }
        set<int> com_neighbors;
        set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
        for (int w : com_neighbors)
        {
            if (find(vertexInfo.changed_edges.begin(), vertexInfo.changed_edges.end(), Edge(v.first, w)) == vertexInfo.changed_edges.end())
            {
                vertexInfo.changed_edges.push_back(Edge(v.first, w));
                vertexInfo.changed_edges.push_back(Edge(w, v.first));
                // if (--kct[w][v.first] < (k - 2))
                if (--kct.at(w).at(v.first) < (k - 2))
                {
                    vertexInfo.sup_less_edges.push_back(make_pair(Edge(w, v.first), 0));
                }
                kct.at(v.first).at(w)--;
            }
        }
    }
    // delete u from kct
    kct.erase(u);
    for (const auto &it : vertexInfo.neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].erase(u);
        }
    }

    // delete the edges whose sup < k-2
    for (int i = 0; i < vertexInfo.sup_less_edges.size(); i++)
    {
        Edge e = vertexInfo.sup_less_edges[i].first;
        int u = e.first;
        int v = e.second;
        vertexInfo.sup_less_edges[i].second = kct.at(u).at(v);
        kct.at(u).erase(v);
        kct.at(v).erase(u);
        // delete isolated vertices
        if (kct[u].empty())
        {
            kct.erase(u);
            del_vertex.insert(u);
        }

        if (kct[v].empty())
        {
            kct.erase(v);
            del_vertex.insert(v);
        }
        if (!kct.count(u) || !kct.count(v))
            continue;

        for (auto nei : kct[u])
        {
            int w = nei.first;
            if (kct.at(v).count(w))
            {
                if (kct[u][w] == (k - 2))
                {
                    vertexInfo.sup_less_edges.push_back(make_pair(Edge(u, w), 0));
                }
                --kct[u][w];
                --kct[w][u];
                if (kct[v][w] == (k - 2))
                {
                    vertexInfo.sup_less_edges.push_back(make_pair(Edge(v, w), 0));
                }
                --kct[v][w];
                --kct[w][v];
            }
        }
    }
    // for (int v : del_vertex)
    // {
    //     kct.erase(v);
    // }
    return move(vertexInfo);
}

void recoverTruss(VertexInfo &vertexInfo, unordered_map<int, unordered_map<int, int>> &kct)
{
    // for (auto &it : vertexInfo.sup_less_edges)
    // {
    //     Edge e = it.first;
    // }
    while (!vertexInfo.sup_less_edges.empty())
    {
        auto &it = vertexInfo.sup_less_edges.back();

        Edge e = it.first;
        int u = e.first;
        int v = e.second;
        if (!kct.count(u))
        {
            kct.emplace(u, unordered_map<int, int>());
        }
        if (!kct.count(v))
        {
            kct.emplace(v, unordered_map<int, int>());
        }
        for (auto nei : kct.at(u))
        {
            int w = nei.first;
            if (kct.at(v).count(w))
            {
                ++kct[u][w];
                ++kct.at(w).at(u);
                ++kct[v][w];
                ++kct.at(w).at(v);
            }
        }
        kct.at(u).emplace(v, it.second);
        kct.at(v).emplace(u, it.second);
        vertexInfo.sup_less_edges.pop_back();
    }
    for (auto &e : vertexInfo.changed_edges)
    {
        kct.at(e.first).at(e.second)++;
    }
    kct.emplace(vertexInfo.vertex, vertexInfo.neighbors);
    for (const auto &it : vertexInfo.neighbors)
    {
        if (kct.count(it.first))
        {
            kct[it.first].emplace(vertexInfo.vertex, it.second);
        }
    }
}
void NaiveEnum_AttTruss(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                        unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            // check connectivity
            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() != M.size())
            {
                return; // 图不连通
            }
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
                record_results(M, score, opt, result);
            }
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "u = " << u << endl;

    set<int> M_temp(M.begin(), M.end());

    // kscore shrink C
    vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att);
    // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    for (int u : removed_vertices)
    {
        if (del_vertex.count(u))
            continue;
        changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
        // del_vertex.insert(u);
    }

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        NaiveEnum_AttTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    for (int v : removed_vertices)
    {
        if (!count(C.begin(), C.end(), v))
            C.push_back(v);
    }

    M.pop_back(); // 恢复 M
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "u = " << u << endl;
    M_temp.erase(u);
    // set<int> bk = {18,6,168,177,205,81,3};
    // set<int> bk2;
    // set<int> bk3 = {513,72,59,219,71,0,125,89,17};
    // set<int> bk4;
    // set_difference(M_temp.begin(), M_temp.end(), bk.begin(), bk.end(), std::inserter(bk2, bk2.begin()));
    // set_difference(M_temp.begin(), M_temp.end(), bk3.begin(), bk3.end(), std::inserter(bk4, bk4.begin()));
    // if(!M_temp.empty()&&bk2.empty()&&bk4.empty()){
    //     int b = 0;
    // }

    // recover kct
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }

    while (!changed_edges.empty())
    {
        auto &it = changed_edges.back();

        recoverTruss(it, kct);
        changed_edges.pop_back();
    }

    changed_edges.clear();
    del_vertex.clear();
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    com.clear();
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        NaiveEnum_AttTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
void NaiveEnum_AttCntTruss(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                           unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
                record_results(M, score, opt, result);
            }
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "u = " << u << endl;

    set<int> M_temp(M.begin(), M.end());

    // kscore shrink C
    vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att);
    // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    for (int u : removed_vertices)
    {
        if (del_vertex.count(u))
            continue;
        changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
        // del_vertex.insert(u);
    }

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        // if (M.empty())
        // {
        //     vector<vector<int>> components = FindConnectedComponents(kct, C);
        //     if (components.size() == 1)
        //     {
        //         NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        //     }
        //     else
        //     {
        //         for (auto &c : components)
        //         {
        //             unordered_map<int, unordered_map<int, int>> local_kct;
        //             for (auto v : c)
        //             {
        //                 local_kct.emplace(v, kct.at(v));
        //             }
        //             NaiveEnum_AttCntTruss(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result);
        //         }
        //     }
        // }
        // else
        // {

        unordered_set<int> visited;
        int startNode = M[0];
        DFS(kct, startNode, visited);
        if (visited.size() == M.size() + C.size())
        {
            NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        }
        else
        {
            vector<int> C_temp;
            int M_count = 0;
            for (int v : visited)
            {
                if (count(M.begin(), M.end(), v))
                    M_count++;
                else
                    C_temp.push_back(v);
            }
            if (M_count == M.size())
            {
                unordered_map<int, unordered_map<int, int>> local_kct;
                for (auto v : visited)
                {
                    local_kct.emplace(v, kct.at(v));
                }
                NaiveEnum_AttCntTruss(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result);
            }
        }
        // }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    for (int v : removed_vertices)
    {
        if (!count(C.begin(), C.end(), v))
            C.push_back(v);
    }

    M.pop_back(); // 恢复 M
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "u = " << u << endl;
    M_temp.erase(u);
    // set<int> bk = {18,6,168,177,205,81,3};
    // set<int> bk2;
    // set<int> bk3 = {513,72,59,219,71,0,125,89,17};
    // set<int> bk4;
    // set_difference(M_temp.begin(), M_temp.end(), bk.begin(), bk.end(), std::inserter(bk2, bk2.begin()));
    // set_difference(M_temp.begin(), M_temp.end(), bk3.begin(), bk3.end(), std::inserter(bk4, bk4.begin()));
    // if(!M_temp.empty()&&bk2.empty()&&bk4.empty()){
    //     int b = 0;
    // }

    // recover kct
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }

    while (!changed_edges.empty())
    {
        auto &it = changed_edges.back();

        recoverTruss(it, kct);
        changed_edges.pop_back();
    }

    changed_edges.clear();
    del_vertex.clear();
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    com.clear();
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
            }
            else
            {
                for (auto &c : components)
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTruss(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result);
                }
            }
        }
        else
        {

            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTruss(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result);
                }
            }
        }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
bool CheckNsize(const vector<int> &M, const vector<int> &C, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, Weight sthd1, Weight sthd2)
{
    // set<int> test = {0, 3, 6, 17, 18, 71, 81, 89, 207, 219};
    // set<int> test2(M.begin(), M.end());
    // set<int> com;
    // set_intersection(test.begin(), test.end(), test2.begin(), test2.end(), std::inserter(com, com.begin()));
    // if (com.size() == test.size())
    // {
    //     int b = 0;
    // }
    int n_M = 0;
    for (auto it : M)
    {
        if (!seq2att[it].first)
        {
            n_M++;
        }
    }
    int n_C = 0;
    for (auto it : C)
    {
        if (!seq2att[it].first)
        {
            n_C++;
        }
    }
    Weight n_C_min = (sthd1 * M.size() - n_M) / (-sthd1 + 1);
    Weight p_C_min = (-sthd2 * M.size() + n_M) / sthd2;

    if (n_C_min > n_C || p_C_min > (C.size() - n_C))
    {
        // cout << "p_C: " << C.size() - n_C << " p_C_min: " << p_C_min << endl;
        return false;
    }
    // if (n_C < n_C_min)
    // {
    //     return false;
    // }
    else
        return true;
}

void NaiveEnum_AttCntTrussNsize(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                                unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         // if (com.size() == result[i].size())
                //         //     result[i] = cand;
                //         // else if (com.size() == cand.size())
                //         // {
                //         // }
                //         // else
                //         //     result.push_back(cand);
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
                record_results(M, score, opt, result);
            }
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    set<int> M_temp(M.begin(), M.end());

    // kscore shrink C
    vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att);
    // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    for (int u : removed_vertices)
    {
        if (del_vertex.count(u))
            continue;
        changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
        // del_vertex.insert(u);
    }

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        // if (M.empty())
        // {
        //     vector<vector<int>> components = FindConnectedComponents(kct);
        //     if (components.size() == 1)
        //     {
        //         NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        //     }
        //     else
        //     {
        //         for (auto &c : components)
        //         {
        //             unordered_map<int, unordered_map<int, int>> local_kct;
        //             for (auto v : c)
        //             {
        //                 local_kct.emplace(v, kct.at(v));
        //             }
        //             NaiveEnum_AttCntTrussNsize(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result);
        //         }
        //     }
        // }
        // else
        // {

        unordered_set<int> visited;
        int startNode = M[0];
        DFS(kct, startNode, visited);
        if (visited.size() == M.size() + C.size())
        {
            if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        }
        else
        {
            vector<int> C_temp;
            int M_count = 0;
            for (int v : visited)
            {
                if (count(M.begin(), M.end(), v))
                    M_count++;
                else
                    C_temp.push_back(v);
            }
            if (M_count == M.size())
            {
                if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTrussNsize(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result);
                }
            }
        }
        // }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    for (int v : removed_vertices)
    {
        if (!count(C.begin(), C.end(), v))
            C.push_back(v);
    }

    M.pop_back(); // 恢复 M
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "u = " << u << endl;
    M_temp.erase(u);
    // set<int> bk = {18,6,168,177,205,81,3};
    // set<int> bk2;
    // set<int> bk3 = {513,72,59,219,71,0,125,89,17};
    // set<int> bk4;
    // set_difference(M_temp.begin(), M_temp.end(), bk.begin(), bk.end(), std::inserter(bk2, bk2.begin()));
    // set_difference(M_temp.begin(), M_temp.end(), bk3.begin(), bk3.end(), std::inserter(bk4, bk4.begin()));
    // if(!M_temp.empty()&&bk2.empty()&&bk4.empty()){
    //     int b = 0;
    // }

    // recover kct
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }

    while (!changed_edges.empty())
    {
        auto &it = changed_edges.back();

        recoverTruss(it, kct);
        changed_edges.pop_back();
    }

    changed_edges.clear();
    del_vertex.clear();
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    com.clear();
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
            }
            else
            {
                for (auto &c : components)
                {

                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTrussNsize(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result);
                }
            }
        }
        else
        {

            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                    NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {
                    if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                    {
                        unordered_map<int, unordered_map<int, int>> local_kct;
                        for (auto v : visited)
                        {
                            local_kct.emplace(v, kct.at(v));
                        }
                        NaiveEnum_AttCntTrussNsize(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result);
                    }
                }
            }
        }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
unordered_map<int, Weight> VertexScore(unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, unordered_map<int, unordered_map<int, int>> &kct, vector<int> &CC)
{
    // vector<int> vert_sort;
    unordered_map<int, Weight> vert_score;
    for (int u : CC)
    {
        Weight min_score = WEIGHT_MAX;
        Weight u_score = seq2att[u].second.first;
        vector<int> u_att = seq2att[u].second.second;
        for (auto it2 : kct[u])
        {
            int v = it2.first;
            Weight v_score = seq2att[v].second.first;
            vector<int> v_att = seq2att[v].second.second;
            vector<int> com_att;
            set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
            Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
            Weight u_v = cal_score(u_score, v_score, u_v_sim);
            min_score = min(min_score, u_v);
        }
        if (min_score != WEIGHT_MAX)
        {
            vert_score[u] = min_score;
        }
        else
        {
            vert_score[u] = WEIGHT_ZERO;
        }
    }
    // vector<pair<int, double>> vec(vert_score.begin(), vert_score.end());

    // // 按第二个元素（double 值）由大到小排序
    // std::sort(vec.begin(), vec.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b)
    //           { return a.second < b.second; });
    // for (auto it : vec)
    // {
    //     vert_sort.push_back(it.first);
    // }
    // return move(vert_sort);
    return vert_score;
}

void NaiveEnum_AttCntTrussNsizeSort(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                                    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k,
                                    Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result, unordered_map<int, Weight> &vert_score)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            Weight score = KScore(M, seq2att, opt);
            if (score > WEIGHT_ZERO)
            {
                // vector<int> cand;
                // for (const auto &it : M)
                // {
                //     cand.push_back(it);
                // }
                // sort(cand.begin(), cand.end());
                // if (score == opt)
                // {
                //     for (int i = 0; i < result.size(); i++)
                //     {
                //         vector<int> com;
                //         set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //         if (com.size() == result[i].size())
                //         {                     //|com|<=|result[i]|
                //             result[i] = cand; // cand>>result[i]
                //             break;
                //         }
                //         else
                //         {
                //             if (com.size() == cand.size()) // result[i]>>cand
                //                 break;
                //             else
                //             {
                //                 if (i == result.size() - 1)
                //                 {
                //                     result.push_back(cand);
                //                 }
                //             }
                //         }
                //     }
                // }
                // else
                // {
                //     result.clear();
                //     result.push_back(cand);
                //     opt = score;
                // }
                record_results(M, score, opt, result);
            }
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    set<int> M_temp(M.begin(), M.end());

    // kscore shrink C
    vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att);
    // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    for (int u : removed_vertices)
    {
        if (del_vertex.count(u))
            continue;
        changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
        // del_vertex.insert(u);
    }

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        // if (M.empty())
        // {
        //     vector<vector<int>> components = FindConnectedComponents(kct);
        //     if (components.size() == 1)
        //     {
        //         NaiveEnum_AttCntTrussNsizeSort(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
        //     }
        //     else
        //     {
        //         for (auto &c : components)
        //         {
        //             unordered_map<int, unordered_map<int, int>> local_kct;
        //             for (auto v : c)
        //             {
        //                 local_kct.emplace(v, kct.at(v));
        //             }
        //             NaiveEnum_AttCntTrussNsizeSort(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
        //         }
        //     }
        // }
        // else
        // {

        unordered_set<int> visited;
        int startNode = M[0];
        DFS(kct, startNode, visited);
        if (visited.size() == M.size() + C.size())
        {
            if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                NaiveEnum_AttCntTrussNsizeSort(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
        }
        else
        {
            vector<int> C_temp;
            int M_count = 0;
            for (int v : visited)
            {
                if (count(M.begin(), M.end(), v))
                    M_count++;
                else
                    C_temp.push_back(v);
            }
            if (M_count == M.size())
            {
                if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTrussNsizeSort(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
                }
            }
        }
        // }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    for (int v : removed_vertices)
    {
        if (!count(C.begin(), C.end(), v))
            C.push_back(v);
    }

    M.pop_back(); // 恢复 M

    M_temp.erase(u);

    // recover kct
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }

    while (!changed_edges.empty())
    {
        auto &it = changed_edges.back();

        recoverTruss(it, kct);
        changed_edges.pop_back();
    }

    changed_edges.clear();
    del_vertex.clear();
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    com.clear();
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_AttCntTrussNsizeSort(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
            }
            else
            {
                for (auto &c : components)
                {

                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTrussNsizeSort(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
                }
            }
        }
        else
        {

            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                    NaiveEnum_AttCntTrussNsizeSort(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {
                    if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                    {
                        unordered_map<int, unordered_map<int, int>> local_kct;
                        for (auto v : visited)
                        {
                            local_kct.emplace(v, kct.at(v));
                        }
                        NaiveEnum_AttCntTrussNsizeSort(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
                    }
                }
            }
        }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
int MaxVertexScore(unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
                   unordered_map<int, unordered_map<int, int>> &kct, vector<int> &M, vector<int> &CC,
                   unordered_map<int, unordered_map<int, Weight>> graph_sim, Weight &thd1, Weight &thd2)
{
    // set<int> test = {17, 207};
    // set<int> test2(CC.begin(), CC.end());
    // set<int> common;
    // set_intersection(test.begin(), test.end(), test2.begin(), test2.end(), std::inserter(common, common.begin()));
    // if (common.size()==2&&common.size() == CC.size())
    //     int a=0;
    int n_M = 0;
    for (auto it : M)
    {
        if (!seq2att[it].first)
        {
            n_M++;
        }
    }
    Weight ns;
    if (M.empty())
    {
        ns = WEIGHT_ZERO;
    }
    else
    {
        ns = Weight(n_M, M.size());
    }

    // vector<int> vert_sort;
    int max_score_vert = -1;
    Weight max_score = WEIGHT_ZERO;
    // double max_score = INT16_MAX;
    // unordered_map<int, double> vert_score;
    int neg = 0;
    if (ns < thd1) // select neg
    {
        neg = -1;
    }
    else if (ns > thd2) // select pos
    {
        neg = 1;
    }
    for (int u : CC)
    {
        if (neg == -1 && seq2att[u].first)
            continue;
        else if (neg == 1 && !seq2att[u].first)
            continue;
        Weight min_score = WEIGHT_MAX;

        for (auto it2 : kct[u])
        {
            int v = it2.first;

            Weight u_v = graph_sim[u][v];
            min_score = min(min_score, u_v);
        }
        if (min_score != WEIGHT_MAX)
        {
            // max_score = max(max_score, min_score);
            if (max_score < min_score)
            {
                max_score = min_score;
                max_score_vert = u;
            }
            else if (max_score == min_score)
            {
                if (kct[max_score_vert].size() < kct[u].size())
                {
                    max_score_vert = u;
                }
            }
        }
    }
    // if(max_score_vert == -1){
    //     int a = 1;
    // }
    return max_score_vert;
}
void NaiveEnum_AttCntTrussNsizeSortNotConst(vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                                            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k,
                                            Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result,
                                            unordered_map<int, unordered_map<int, Weight>> &graph_sim)
{
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            Weight score = KScore(M, seq2att, opt);
            if (score > 0)
            {
                //     vector<int> cand;
                //     for (const auto &it : M)
                //     {
                //         cand.push_back(it);
                //     }
                //     sort(cand.begin(), cand.end());
                //     if (score == opt)
                //     {
                //         if (result.size() == 0)
                //         {
                //             result.push_back(cand);
                //             result_score.push_back(score);
                //         }
                //         else
                //         {
                //             for (int i = 0; i < result.size(); i++)
                //             {
                //                 vector<int> com;
                //                 set_intersection(cand.begin(), cand.end(), result[i].begin(), result[i].end(), back_inserter(com));
                //                 if (com.size() == result[i].size())
                //                 {                     //|com|<=|result[i]|
                //                     result[i] = cand; // cand>>result[i]
                //                     result_score[i] = score;
                //                     // cout << "cand score: " << score << endl;
                //                     break;
                //                 }
                //                 else
                //                 {
                //                     if (com.size() == cand.size()) // result[i]>>cand
                //                         break;
                //                     else
                //                     {
                //                         if (i == result.size() - 1)
                //                         {
                //                             result.push_back(cand);
                //                             result_score.push_back(score);
                //                             // cout << "cand score: " << score << endl;
                //                         }
                //                     }
                //                 }
                //             }
                //         }
                //     }
                //     else
                //     {
                //         result.clear();
                //         result.push_back(cand);
                //         result_score.clear();
                //         result_score.push_back(score);
                //         // cout << "cand score: " << score << endl;
                //         opt = score;
                //     }
                record_results(M, score, opt, result);
            }
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    // int u = C.back();
    // C.pop_back();
    int u = MaxVertexScore(seq2att, kct, M, C, graph_sim, sthd1, sthd2); // 可以优化，这里每个搜索点都计算一次，但很多分支的kct不变，导致重复计算
    C.erase(std::find(C.begin(), C.end(), u));

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    set<int> M_temp(M.begin(), M.end());

    // kscore shrink C
    vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att);
    // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    for (int u : removed_vertices)
    {
        if (del_vertex.count(u))
            continue;
        changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
        // del_vertex.insert(u);
    }

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        // if (M.empty())
        // {
        //     vector<vector<int>> components = FindConnectedComponents(kct);
        //     if (components.size() == 1)
        //     {
        //         NaiveEnum_AttCntTrussNsizeSortNotConst(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, result_score);
        //     }
        //     else
        //     {
        //         for (auto &c : components)
        //         {
        //             unordered_map<int, unordered_map<int, int>> local_kct;
        //             for (auto v : c)
        //             {
        //                 local_kct.emplace(v, kct.at(v));
        //             }
        //             NaiveEnum_AttCntTrussNsizeSortNotConst(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, result_score);
        //         }
        //     }
        // }
        // else
        // {

        unordered_set<int> visited;
        int startNode = M[0];
        DFS(kct, startNode, visited);
        if (visited.size() == M.size() + C.size())
        {
            if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                NaiveEnum_AttCntTrussNsizeSortNotConst(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
        }
        else
        {
            vector<int> C_temp;
            int M_count = 0;
            for (int v : visited)
            {
                if (count(M.begin(), M.end(), v))
                    M_count++;
                else
                    C_temp.push_back(v);
            }
            if (M_count == M.size())
            {
                if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTrussNsizeSortNotConst(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
                }
            }
        }
        // }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    for (int v : removed_vertices)
    {
        if (!count(C.begin(), C.end(), v))
            C.push_back(v);
    }

    M.pop_back(); // 恢复 M

    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;

    M_temp.erase(u);

    // recover kct
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }

    while (!changed_edges.empty())
    {
        auto &it = changed_edges.back();

        recoverTruss(it, kct);
        changed_edges.pop_back();
    }

    changed_edges.clear();
    del_vertex.clear();
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    com.clear();
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_AttCntTrussNsizeSortNotConst(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
            }
            else
            {
                for (auto &c : components)
                {

                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTrussNsizeSortNotConst(M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
                }
            }
        }
        else
        {
            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                    NaiveEnum_AttCntTrussNsizeSortNotConst(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {
                    if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                    {
                        unordered_map<int, unordered_map<int, int>> local_kct;
                        for (auto v : visited)
                        {
                            local_kct.emplace(v, kct.at(v));
                        }
                        NaiveEnum_AttCntTrussNsizeSortNotConst(M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
                    }
                }
            }
        }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
void NaiveEnum_AttCntTrussNsizeSortNotConst(Weight cur_score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                                            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k,
                                            Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result,
                                            unordered_map<int, unordered_map<int, Weight>> &graph_sim, chrono::high_resolution_clock::time_point startTime)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - startTime;
    if (duration.count() > max_time)
    { // 超过1秒，终止递归
        // std::cout << "Timeout: Execution exceeded " << max_time << " second." << std::endl;
        result.clear();
        return;
    }
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            record_results(M, cur_score, opt, result);
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();
    // if (cur_score < opt)
    //     return;
    // int u = MaxVertexScore(seq2att, kct, M, C, graph_sim, sthd1, sthd2);
    // C.erase(std::find(C.begin(), C.end(), u));
    Weight min_score = WEIGHT_MAX;
    // Weight u_score = seq2att[u].second.first;
    // vector<int> u_att = seq2att[u].second.second;
    for (int v : M)
    {
        // Weight v_score = seq2att[v].second.first;
        // vector<int> v_att = seq2att[v].second.second;
        // vector<int> com_att;
        // set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
        // Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
        // Weight u_v = cal_score(u_score, v_score, u_v_sim);
        // if (u_v < min_score)
        // {
        //     min_score = u_v;
        // }
        Weight u_v = graph_sim[u][v];
        if (u_v < min_score)
        {
            min_score = u_v;
        }
    }
    set<int> M_temp(M.begin(), M.end());
    if (min_score >= opt)
    { // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
        M.push_back(u);
        M_temp.insert(u);
        // for (auto &it : M)
        // {
        //     cout << it << " ";
        // }
        // cout << ", ";
        // for (auto &it : C)
        // {
        //     cout << it << " ";
        // }
        // cout << endl;

        // kscore shrink C
        vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att, graph_sim);
        // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
        vector<VertexInfo> changed_edges;
        set<int> del_vertex;
        for (int u : removed_vertices)
        {
            if (del_vertex.count(u))
                continue;
            changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
            // del_vertex.insert(u);
        }

        set<int> com;
        set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
        if (com.size() == 0)
        {
            for (int v : del_vertex)
            {
                auto it = std::find(C.begin(), C.end(), v);
                if (it != C.end())
                {
                    // 如果找到，则删除
                    C.erase(it);
                }
            }
            //////////////////////////////////////

            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                    NaiveEnum_AttCntTrussNsizeSortNotConst(min(cur_score, min_score), M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {
                    if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                    {
                        unordered_map<int, unordered_map<int, int>> local_kct;
                        for (auto v : visited)
                        {
                            local_kct.emplace(v, kct.at(v));
                        }
                        NaiveEnum_AttCntTrussNsizeSortNotConst(min(cur_score, min_score), M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
                    }
                }
            }
            // }
            //////////////////////////////////////
            // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
            for (int v : del_vertex)
            {
                C.push_back(v);
            }
        }
        for (int v : removed_vertices)
        {
            if (!count(C.begin(), C.end(), v))
                C.push_back(v);
        }
        // recover kct
        for (int v : del_vertex)
        {
            kct.emplace(v, unordered_map<int, int>());
        }

        while (!changed_edges.empty())
        {
            auto &it = changed_edges.back();

            recoverTruss(it, kct);
            changed_edges.pop_back();
        }
        M.pop_back(); // 恢复 M
    }
    if (cur_score < opt)
        return;
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "0" << endl;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_AttCntTrussNsizeSortNotConst(cur_score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
            }
            else
            {
                for (auto &c : components)
                {

                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttCntTrussNsizeSortNotConst(cur_score, M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
                }
            }
        }
        else
        {
            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                    NaiveEnum_AttCntTrussNsizeSortNotConst(cur_score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {
                    if (CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                    {
                        unordered_map<int, unordered_map<int, int>> local_kct;
                        for (auto v : visited)
                        {
                            local_kct.emplace(v, kct.at(v));
                        }
                        NaiveEnum_AttCntTrussNsizeSortNotConst(cur_score, M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
                    }
                }
            }
        }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
void NaiveEnum_AttTrussCnnt(Weight cur_score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k,
                            Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result,
                            unordered_map<int, unordered_map<int, Weight>> &graph_sim, chrono::high_resolution_clock::time_point startTime)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - startTime;
    if (duration.count() > max_time)
    { // 超过1秒，终止递归
        // std::cout << "Timeout: Execution exceeded " << max_time << " second." << std::endl;
        result.clear();
        return;
    }
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            record_results(M, cur_score, opt, result);
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();
    // if (cur_score < opt)
    //     return;
    // int u = MaxVertexScore(seq2att, kct, M, C, graph_sim, sthd1, sthd2);
    // C.erase(std::find(C.begin(), C.end(), u));
    Weight min_score = WEIGHT_MAX;
    // Weight u_score = seq2att[u].second.first;
    // vector<int> u_att = seq2att[u].second.second;
    for (int v : M)
    {
        // Weight v_score = seq2att[v].second.first;
        // vector<int> v_att = seq2att[v].second.second;
        // vector<int> com_att;
        // set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
        // Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
        // Weight u_v = cal_score(u_score, v_score, u_v_sim);
        // if (u_v < min_score)
        // {
        //     min_score = u_v;
        // }
        Weight u_v = graph_sim[u][v];
        if (u_v < min_score)
        {
            min_score = u_v;
        }
    }
    set<int> M_temp(M.begin(), M.end());
    if (min_score >= opt)
    { // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
        M.push_back(u);
        M_temp.insert(u);
        // for (auto &it : M)
        // {
        //     cout << it << " ";
        // }
        // cout << ", ";
        // for (auto &it : C)
        // {
        //     cout << it << " ";
        // }
        // cout << endl;

        // kscore shrink C
        vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att, graph_sim);
        // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
        vector<VertexInfo> changed_edges;
        set<int> del_vertex;
        for (int u : removed_vertices)
        {
            if (del_vertex.count(u))
                continue;
            changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
            // del_vertex.insert(u);
        }

        set<int> com;
        set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
        if (com.size() == 0)
        {
            for (int v : del_vertex)
            {
                auto it = std::find(C.begin(), C.end(), v);
                if (it != C.end())
                {
                    // 如果找到，则删除
                    C.erase(it);
                }
            }
            //////////////////////////////////////

            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                NaiveEnum_AttTrussCnnt(min(cur_score, min_score), M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {

                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttTrussCnnt(min(cur_score, min_score), M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
                }
            }
            // }
            //////////////////////////////////////
            // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
            for (int v : del_vertex)
            {
                C.push_back(v);
            }
        }
        for (int v : removed_vertices)
        {
            if (!count(C.begin(), C.end(), v))
                C.push_back(v);
        }
        // recover kct
        for (int v : del_vertex)
        {
            kct.emplace(v, unordered_map<int, int>());
        }

        while (!changed_edges.empty())
        {
            auto &it = changed_edges.back();

            recoverTruss(it, kct);
            changed_edges.pop_back();
        }
        M.pop_back(); // 恢复 M
    }
    if (cur_score < opt)
        return;
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "0" << endl;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_AttTrussCnnt(cur_score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
            }
            else
            {
                for (auto &c : components)
                {

                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttTrussCnnt(cur_score, M, c, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
                }
            }
        }
        else
        {
            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size())
            {
                NaiveEnum_AttTrussCnnt(cur_score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
            }
            else
            {
                vector<int> C_temp;
                int M_count = 0;
                for (int v : visited)
                {
                    if (count(M.begin(), M.end(), v))
                        M_count++;
                    else
                        C_temp.push_back(v);
                }
                if (M_count == M.size())
                {

                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_AttTrussCnnt(cur_score, M, C_temp, local_kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
                }
            }
        }
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
void NaiveEnum_AttTruss(Weight cur_score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                        unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k,
                        Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result,
                        unordered_map<int, unordered_map<int, Weight>> &graph_sim, chrono::high_resolution_clock::time_point startTime)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - startTime;
    if (duration.count() > max_time)
    { // 超过1秒，终止递归
        // std::cout << "Timeout: Execution exceeded " << max_time << " second." << std::endl;
        result.clear();
        return;
    }
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        if (M.empty())
            return;
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            // check connectivity
            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() != M.size())
            {
                return; // 图不连通
            }
            record_results(M, cur_score, opt, result);
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();
    // if (cur_score < opt)
    //     return;
    // int u = MaxVertexScore(seq2att, kct, M, C, graph_sim, sthd1, sthd2);
    // C.erase(std::find(C.begin(), C.end(), u));
    Weight min_score = WEIGHT_MAX;
    // Weight u_score = seq2att[u].second.first;
    // vector<int> u_att = seq2att[u].second.second;
    for (int v : M)
    {
        // Weight v_score = seq2att[v].second.first;
        // vector<int> v_att = seq2att[v].second.second;
        // vector<int> com_att;
        // set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
        // Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
        // Weight u_v = cal_score(u_score, v_score, u_v_sim);
        // if (u_v < min_score)
        // {
        //     min_score = u_v;
        // }
        Weight u_v = graph_sim[u][v];
        if (u_v < min_score)
        {
            min_score = u_v;
        }
    }
    set<int> M_temp(M.begin(), M.end());
    if (min_score >= opt)
    { // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
        M.push_back(u);
        M_temp.insert(u);
        // for (auto &it : M)
        // {
        //     cout << it << " ";
        // }
        // cout << ", ";
        // for (auto &it : C)
        // {
        //     cout << it << " ";
        // }
        // cout << endl;

        // kscore shrink C
        vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att, graph_sim);
        // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
        vector<VertexInfo> changed_edges;
        set<int> del_vertex;
        for (int u : removed_vertices)
        {
            if (del_vertex.count(u))
                continue;
            changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
            // del_vertex.insert(u);
        }

        set<int> com;
        set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
        if (com.size() == 0)
        {
            for (int v : del_vertex)
            {
                auto it = std::find(C.begin(), C.end(), v);
                if (it != C.end())
                {
                    // 如果找到，则删除
                    C.erase(it);
                }
            }
            //////////////////////////////////////

            NaiveEnum_AttTruss(min(cur_score, min_score), M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);

            //////////////////////////////////////
            // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
            for (int v : del_vertex)
            {
                C.push_back(v);
            }
        }
        for (int v : removed_vertices)
        {
            if (!count(C.begin(), C.end(), v))
                C.push_back(v);
        }
        // recover kct
        for (int v : del_vertex)
        {
            kct.emplace(v, unordered_map<int, int>());
        }

        while (!changed_edges.empty())
        {
            auto &it = changed_edges.back();

            recoverTruss(it, kct);
            changed_edges.pop_back();
        }
        M.pop_back(); // 恢复 M
    }
    if (cur_score < opt)
        return;
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "0" << endl;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////

        NaiveEnum_AttTruss(cur_score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
void NaiveEnum_AttTrussSize(Weight cur_score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k,
                            Weight &opt, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result,
                            unordered_map<int, unordered_map<int, Weight>> &graph_sim, chrono::high_resolution_clock::time_point startTime)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - startTime;
    if (duration.count() > max_time)
    { // 超过1秒，终止递归
        // std::cout << "Timeout: Execution exceeded " << max_time << " second." << std::endl;
        result.clear();
        return;
    }
    // result stores the result whose score is larger than opt
    if (C.empty())
    {
        if (M.empty())
            return;
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            // check connectivity
            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() != M.size())
            {
                return; // 图不连通
            }
            record_results(M, cur_score, opt, result);
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();
    // if (cur_score < opt)
    //     return;
    // int u = MaxVertexScore(seq2att, kct, M, C, graph_sim, sthd1, sthd2);
    // C.erase(std::find(C.begin(), C.end(), u));
    Weight min_score = WEIGHT_MAX;
    // Weight u_score = seq2att[u].second.first;
    // vector<int> u_att = seq2att[u].second.second;
    for (int v : M)
    {
        // Weight v_score = seq2att[v].second.first;
        // vector<int> v_att = seq2att[v].second.second;
        // vector<int> com_att;
        // set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
        // Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
        // Weight u_v = cal_score(u_score, v_score, u_v_sim);
        // if (u_v < min_score)
        // {
        //     min_score = u_v;
        // }
        Weight u_v = graph_sim[u][v];
        if (u_v < min_score)
        {
            min_score = u_v;
        }
    }
    set<int> M_temp(M.begin(), M.end());
    if (min_score >= opt)
    { // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
        M.push_back(u);
        M_temp.insert(u);
        // for (auto &it : M)
        // {
        //     cout << it << " ";
        // }
        // cout << ", ";
        // for (auto &it : C)
        // {
        //     cout << it << " ";
        // }
        // cout << endl;

        // kscore shrink C
        vector<int> removed_vertices = Shrink_C_Att(u, opt, C, seq2att, graph_sim);
        // unordered_map<int, unordered_map<int, int>> kct_temp = kct;
        vector<VertexInfo> changed_edges;
        set<int> del_vertex;
        for (int u : removed_vertices)
        {
            if (del_vertex.count(u))
                continue;
            changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
            // del_vertex.insert(u);
        }

        set<int> com;
        set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
        if (com.size() == 0)
        {
            for (int v : del_vertex)
            {
                auto it = std::find(C.begin(), C.end(), v);
                if (it != C.end())
                {
                    // 如果找到，则删除
                    C.erase(it);
                }
            }
            //////////////////////////////////////
            if (CheckNsize(M, C, seq2att, sthd1, sthd2))
                NaiveEnum_AttTrussSize(min(cur_score, min_score), M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);

            //////////////////////////////////////
            // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
            for (int v : del_vertex)
            {
                C.push_back(v);
            }
        }
        for (int v : removed_vertices)
        {
            if (!count(C.begin(), C.end(), v))
                C.push_back(v);
        }
        // recover kct
        for (int v : del_vertex)
        {
            kct.emplace(v, unordered_map<int, int>());
        }

        while (!changed_edges.empty())
        {
            auto &it = changed_edges.back();

            recoverTruss(it, kct);
            changed_edges.pop_back();
        }
        M.pop_back(); // 恢复 M
    }
    if (cur_score < opt)
        return;
    // for (auto &it : M)
    // {
    //     cout << it << " ";
    // }
    // cout << ", ";
    // for (auto &it : C)
    // {
    //     cout << it << " ";
    // }
    // cout << endl;
    // cout << "0" << endl;
    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    // kct = kct_temp;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));

    set<int> com;
    set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
    if (com.size() == 0)
    {
        for (int v : del_vertex)
        {
            auto it = std::find(C.begin(), C.end(), v);
            if (it != C.end())
            {
                // 如果找到，则删除
                C.erase(it);
            }
        }
        //////////////////////////////////////
        if (CheckNsize(M, C, seq2att, sthd1, sthd2))
            NaiveEnum_AttTrussSize(cur_score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
        //////////////////////////////////////
        // NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd, result);
        for (int v : del_vertex)
        {
            C.push_back(v);
        }
    }
    // 恢复 C
    C.push_back(u);
    // kct = kct_temp;
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
void GetKTrussInc(unordered_map<int, unordered_map<int, int>> &graph, vector<Edge> &sup_less_edges, vector<Edge> &added_graph, int k)
{
    // recalculate the added_graph(contain the added edges and sup_less_edges)
    for (auto &it : added_graph)
    {

        update_sup(graph, it);
    }
    for (auto &it : sup_less_edges)
    {
        update_sup(graph, it);
    }

    // unordered_map<int, unordered_map<int, int>> ktruss;
    sup_less_edges.clear();
    vector<Edge> temp_sup_less_edges;
    for (auto &it : graph)
    {
        int u = it.first;
        for (auto &it1 : it.second)
        {
            int v = it1.first;
            if (v > u)
                continue;
            int sup = it1.second;
            if (sup < k - 2)
            {
                temp_sup_less_edges.push_back(pair<int, int>(u, v));
            }
        }
    }
    MaintainKTruss_DeleteEdge(graph, sup_less_edges, temp_sup_less_edges, k);
}

void DFS(const unordered_map<int, unordered_map<int, int>> &graph, int start, unordered_set<int> &visited, const unordered_set<int> &neighbors)
{
    stack<int> stack;
    stack.push(start);
    while (!stack.empty())
    {
        int node = stack.top();
        stack.pop();
        if (visited.find(node) == visited.end())
        {
            visited.insert(node);
            for (const auto &neighbor : graph.at(node))
            {
                if (neighbors.find(neighbor.first) != neighbors.end())
                {
                    stack.push(neighbor.first);
                }
            }
        }
    }
}

// 检查子图是否连通的函数
bool isSubgraphConnected(const unordered_map<int, unordered_map<int, int>> &graph, const unordered_set<int> &neighbors)
{
    if (neighbors.empty())
    {
        return true; // 空子图被认为是连通的
    }

    unordered_set<int> visited;
    // 选择 neighbors 中的一个节点作为起始节点
    int startNode = *neighbors.begin();

    // 从起始节点开始进行 DFS 遍历
    DFS(graph, startNode, visited, neighbors);

    // 检查 neighbors 中的所有节点是否都被访问过
    if (visited.size() == neighbors.size())
        return true; // 所有节点都被访问过，子图连通
    else
        return false;
}
// double calcul_weight(unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int v,
//                      unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     double weight = 0.0;
//     for (auto &it2 : graph[v])
//     {
//         int u = it2.first;
//         // if (count(cc.begin(), cc.end(), u))
//         // {
//         // weight += 1.0 / (1.0 + 10 * it2.second);
//         Weight w = sim_graph[u][v] + 1;
//         weight = weight + (log(w.numerator) - log(w.denominator));
//         // weight = weight + (log(1 + it2.second.numerator) - log(1 + it2.second.denominator));
//         // }
//     }
//     // if (sim_neighbors > 0)
//     // {
//     //     weight = -(weight / sim_neighbors); // 计算平均权重
//     //     // weight = exp(weight / sim_neighbors);
//     // }
//     // else
//     // {
//     //     weight = 0.0; // 如果没有有效邻居，设置 weight 为 0
//     // }
//     return weight / graph[v].size();
// }
double calcul_weight(unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int v,
                     unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    double weight = 0.0;
    for (auto &it2 : sim_graph[v])
    {
        int u = it2.first;
        if (count(cc.begin(), cc.end(), u))
        {
            // weight += 1.0 / (1.0 + 10 * it2.second);
            // Weight w = sim_graph[u][v] + 1;
            // weight = weight + (log(w.numerator) - log(w.denominator));
            // Weight w = sim_graph[u][v] + 1;
            Weight w = sim_graph[u][v];
            weight += double(w.numerator) / w.denominator;
            // weight = weight + (log(1 + it2.second.numerator) - log(1 + it2.second.denominator));
        }
    }
    // if (sim_neighbors > 0)
    // {
    //     weight = -(weight / sim_neighbors); // 计算平均权重
    //     // weight = exp(weight / sim_neighbors);
    // }
    // else
    // {
    //     weight = 0.0; // 如果没有有效邻居，设置 weight 为 0
    // }
    // return log(weight) - log(cc.size());
    return weight / (cc.size() - 1);
}
// Weight calcul_weight(unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int v,
//                      unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     Weight weight(0, 1);
//     int sim_neighbors = 0;
//     for (auto &it2 : sim_graph[v])
//     {
//         int u = it2.first;
//         if (count(cc.begin(), cc.end(), u))
//         {
//             sim_neighbors++;
//             // weight += 1.0 / (1.0 + 10 * it2.second);
//             weight = weight + it2.second;
//         }
//     }
//     if (sim_neighbors > 0)
//     {
//         weight = weight / sim_neighbors; // 计算平均权重
//         // weight = exp(weight / sim_neighbors);
//     }
//     // else
//     // {
//     //     weight = 0.0; // 如果没有有效邻居，设置 weight 为 0
//     // }
//     return weight;
// }
// double calcul_weight(unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int v,
//                      unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     double weight = 0.0;
//     int sim_neighbors = 0;
//     for (auto &it2 : sim_graph[v])
//     {
//         int u = it2.first;
//         if (count(cc.begin(), cc.end(), u))
//         {
//             sim_neighbors++;
//             // weight += 1.0 / (1.0 + 10 * it2.second);
//             weight = weight + (log(it2.second.numerator) - log(it2.second.denominator));
//         }
//     }
//     if (sim_neighbors > 0)
//     {
//         // weight = weight / sim_neighbors; // 计算平均权重
//         weight = exp(weight / sim_neighbors);
//     }
//     // else
//     // {
//     //     weight = 0.0; // 如果没有有效邻居，设置 weight 为 0
//     // }
//     return weight;
// }
// double calcul_weight(unordered_map<int, unordered_map<int, int>> &graph, int v, unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     double weight = 1;
//     int sim_neighbors = 0;
//     for (auto &it2 : sim_graph[v])
//     {
//         int u = it2.first;
//         if (graph.count(u))
//         {
//             sim_neighbors++;
//             // weight += 1.0 / (1.0 + 10 * it2.second);
//             weight = weight * log(it2.second.numerator) / log(it2.second.denominator);
//         }
//     }
//     if (sim_neighbors > 0)
//     {
//         weight = weight / sim_neighbors; // 计算平均权重
//     }
//     // else
//     // {
//     //     weight = 0.0; // 如果没有有效邻居，设置 weight 为 0
//     // }
//     return weight;
// }
bool compareByValue(const pair<int, double> &a, const pair<int, double> &b)
{
    return a.second < b.second;
}
// bool compareByValue(const pair<int, Weight> &a, const pair<int, Weight> &b)
// {
//     return a.second < b.second;
// }
void isEdgedVertex(unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int v,
                   unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
                   unordered_map<int, double> &pos_edged_vertices, unordered_map<int, double> &neg_edged_vertices,
                   unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    unordered_set<int> neighbors;
    for (auto &it2 : graph[v])
    {
        int dst = it2.first;
        neighbors.insert(dst);
    }
    if (isSubgraphConnected(graph, neighbors))
    {
        double weight = calcul_weight(graph, cc, v, sim_graph);
        // Weight weight = calcul_weight(graph, cc, v, sim_graph);
        if (seq2att[v].first)
        {
            pos_edged_vertices.emplace(v, weight);
        }
        else
        {
            // neg_edged_vertices.emplace(v, sup);
            neg_edged_vertices.emplace(v, weight);
        }
    }
}
// judge whether the vertex can be deleted.
// bool Check_Edge_vertex(int u, unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int k,
//                        unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, vector<VertexInfo> &Removed_vertices,
//                        unordered_map<int, double> &pos_edged_vertices, unordered_map<int, double> &neg_edged_vertices,
//                        unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     bool flag = true;
//     vector<Edge> changed_edges2;
//     unordered_map<int, int> neighbors; // the new edged vertex must be the neighbor of the removed vertex.
//     // visit all opposite edges of u. there are two case: 1. all opposite edges are not deleted. 2. a opposite edge is deleted.
//     for (auto it = graph[u].begin(); it != graph[u].end();)
//     {
//         int v = it->first;

//         for (auto &it1 : graph[v])
//         {
//             int w = it1.first;
//             // int new_sup = -1;
//             // update support
//             // 需要证明如果边缘点候选存在临边是被删除检查点影响的边，那么检查点的删除不影响这个点是否还是边缘点。具体的，他们的关系是共同构成一个团
//             if (graph[u].count(w))
//             {
//                 if (graph[v][w] == k - 2)
//                 {
//                     flag = false;
//                     break;
//                 }
//                 --graph[v][w];
//                 --graph[w][v];
//                 changed_edges2.push_back(pair<int, int>(v, w));
//             }
//         }
//         if (!flag)
//             break;
//         neighbors.emplace(v, it->second);
//         it = graph[u].erase(it);
//     }
//     if (flag)
//     {
//         cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
//         graph.erase(u);
//         vector<int> new_vertex;
//         for (auto &v_sup : neighbors)
//         {
//             int v = v_sup.first;
//             graph[v].erase(u); // remove u from v's neighbors
//             // check whether has new edged vertex
//             if (pos_edged_vertices.count(v) || neg_edged_vertices.count(v))
//             {
//                 Weight a = sim_graph[u][v];
//                 auto it = pos_edged_vertices.find(v);
//                 if (it != pos_edged_vertices.end())
//                 {

//                     it->second -= log(1 + a.numerator) - log(1 + a.denominator);
//                     // update_weight(graph, u, v, sim_graph);
//                 }
//                 else
//                 {
//                     // neg_edged_vertices[v] = calcul_weight(graph, cc, v, sim_graph);
//                     neg_edged_vertices[v] -= log(1 + a.numerator) - log(1 + a.denominator);
//                 }
//                 // continue;
//             }
//             else
//             {
//                 new_vertex.push_back(v);
//                 // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//             }
//         }
//         for (auto &v : new_vertex)
//         {
//             isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//         }
//         VertexInfo removeV;
//         removeV.vertex = u;
//         removeV.neighbors = move(neighbors);
//         removeV.changed_edges = move(changed_edges2);
//         Removed_vertices.push_back(removeV);

//         // update ns

//         return true;
//     }
//     else
//     {
//         // recovery
//         for (auto &v_sup : neighbors)
//         {
//             graph[u].emplace(v_sup.first, v_sup.second);
//         }
//         for (auto &it : changed_edges2)
//         {
//             int v = it.first;
//             int w = it.second;
//             ++graph[v][w];
//             ++graph[w][v];
//         }
//         return false;
//     }
// }
// 5-20
bool Check_Edge_vertex(int u, unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int k,
                       unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, vector<VertexInfo> &Removed_vertices,
                       unordered_map<int, double> &pos_edged_vertices, unordered_map<int, double> &neg_edged_vertices,
                       unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    bool flag = true;
    vector<Edge> changed_edges2;
    unordered_map<int, int> neighbors; // the new edged vertex must be the neighbor of the removed vertex.
    // visit all opposite edges of u. there are two case: 1. all opposite edges are not deleted. 2. a opposite edge is deleted.
    for (auto it = graph[u].begin(); it != graph[u].end();)
    {
        int v = it->first;

        for (auto &it1 : graph[v])
        {
            int w = it1.first;
            // int new_sup = -1;
            // update support
            // 需要证明如果边缘点候选存在临边是被删除检查点影响的边，那么检查点的删除不影响这个点是否还是边缘点。具体的，他们的关系是共同构成一个团
            if (graph[u].count(w))
            {
                if (graph[v][w] == k - 2)
                {
                    flag = false;
                    break;
                }
                --graph[v][w];
                --graph[w][v];
                changed_edges2.push_back(pair<int, int>(v, w));
            }
        }
        if (!flag)
            break;
        neighbors.emplace(v, it->second);
        it = graph[u].erase(it);
    }
    if (flag)
    {
        for (auto &v_w : pos_edged_vertices)
        {

            int v = v_w.first;
            if (u == v)
                continue;
            Weight a = sim_graph[u][v] + 1;
            v_w.second -= log(a.numerator) - log(a.denominator);
        }
        for (auto &v_w : neg_edged_vertices)
        {
            int v = v_w.first;
            if (u == v)
                continue;
            Weight a = sim_graph[u][v] + 1;
            v_w.second -= log(a.numerator) - log(a.denominator);
        }
        cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
        graph.erase(u);
        vector<int> new_vertex;
        for (auto &v_sup : neighbors)
        {
            int v = v_sup.first;
            graph[v].erase(u); // remove u from v's neighbors

            // check whether has new edged vertex
            if (!pos_edged_vertices.count(v) && !neg_edged_vertices.count(v))
            {
                new_vertex.push_back(v);
                // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
            }
        }
        for (auto &v : new_vertex)
        {
            isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
        }
        VertexInfo removeV;
        removeV.vertex = u;
        removeV.neighbors = move(neighbors);
        removeV.changed_edges = move(changed_edges2);
        Removed_vertices.push_back(removeV);

        // update ns

        return true;
    }
    else
    {
        // recovery
        for (auto &v_sup : neighbors)
        {
            graph[u].emplace(v_sup.first, v_sup.second);
        }
        for (auto &it : changed_edges2)
        {
            int v = it.first;
            int w = it.second;
            ++graph[v][w];
            ++graph[w][v];
        }
        return false;
    }
}
void Deleted_vertex(int u, unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int k,
                    unordered_set<int> &deleted_vertices,
                    unordered_map<int, unordered_map<int, int>> &deleted_edges)
{
    // get 1-hop subgraph
    unordered_map<int, unordered_map<int, int>> subgraph;
    const auto &neighbors = graph[u];
    for (const auto &[v, weight_v] : neighbors)
    {
        subgraph[v] = {};
        // subgraph[u][v] = weight_v; // 添加 u -> v 边
        for (const auto &[w, weight_w] : graph.at(v))
        {
            // 只保留邻居之间的连接
            if (neighbors.count(w))
            {
                subgraph[v][w] = weight_w;
            }
        }
    }
    unordered_set<int> deleted_vertices_cand;
    for (auto &u_v_w : subgraph)
    {
        int u = u_v_w.first;
        if (subgraph[u].size() == graph[u].size() - 1)
        {
            deleted_vertices_cand.insert(u);
        }
    }
    // unordered_map<int, unordered_map<int, int>> deleted_edges;
    for (auto it = neighbors.begin(); it != neighbors.end(); it++)
    {
        int v = it->first;

        for (auto &it1 : subgraph[v])
        {

            int w = it1.first;
            if (v > w)
                continue;
            // int new_sup = -1;
            // update support
            // 需要证明如果边缘点候选存在临边是被删除检查点影响的边，那么检查点的删除不影响这个点是否还是边缘点。具体的，他们的关系是共同构成一个团
            int sup = subgraph[v][w];
            if (sup == k - 2)
            {
                if (!deleted_vertices_cand.count(w) && !deleted_vertices_cand.count(v))
                {
                    return;
                }
                if (!deleted_edges.count(v))
                {
                    deleted_edges.emplace(v, unordered_map<int, int>());
                    // deleted_edges.emplace(v, unordered_set<int>());
                }
                deleted_edges[v].emplace(w, sup);
                // deleted_edges[v].insert(w);
                if (!deleted_edges.count(w))
                {
                    deleted_edges.emplace(w, unordered_map<int, int>());
                    // deleted_edges.emplace(w, unordered_set<int>());
                }
                deleted_edges[w].emplace(v, sup);
                // deleted_edges[w].insert(v);
            }
            // else
            // {
            //     changed_edges.push_back(Edge(v, w));
            // }
        }
    }
    // unordered_set<int> deleted_vertices;
    for (auto &u_v : deleted_edges)
    {
        int u = u_v.first;
        // if (subgraph[u].size() != graph[u].size() - 1)
        // {
        //     deleted_vertices.clear();
        //     // deleted_edges.clear();
        //     return;
        // }
        if (deleted_vertices_cand.count(u) && deleted_edges[u].size() == subgraph[u].size())
        {
            deleted_vertices.insert(u);
        }
    }

    for (auto &u : deleted_vertices)
    {

        for (auto &v_s : subgraph[u])
        {
            int v = v_s.first;
            subgraph[v].erase(u);
        }
        subgraph.erase(u);
    }
    if (subgraph.empty())
    {
        deleted_vertices.clear();
        deleted_edges.clear();
        return;
    }
    unordered_set<int> visited;
    int w = subgraph.begin()->first;
    DFS(subgraph, w, visited);
    if (visited.size() == subgraph.size())
    {
        deleted_vertices.insert(u);
        // deleted_edges.emplace(u, unordered_map<int, int>());

        // for (auto &v_sup : graph[u])
        // {
        //     int v = v_sup.first;
        //     int sup = v_sup.second;
        //     deleted_edges[u].emplace(v, sup);
        //     if (!deleted_edges.count(v))
        //     {
        //         deleted_edges.emplace(v, unordered_map<int, int>());
        //         // deleted_edges.emplace(v, unordered_set<int>());
        //     }
        //     deleted_edges[v].emplace(w, sup);
        //     // deleted_edges[v].insert(w);
        // }
    }
    else
    {
        deleted_vertices.clear();
        deleted_edges.clear();
        // changed_edges.clear();
    }
    // if (visited.size() == subgraph.size())
    // {
    //     // deleted_vertices.insert(u);
    //     // return move(deleted_vertices);
    // }
    // else
    // {
    //     return {};
    // }
}

bool First_Pos_Edged_Vertex(unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int k, Weight &ns,
                            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, vector<VertexInfo> &Removed_vertices,
                            unordered_map<int, double> &pos_edged_vertices, unordered_map<int, double> &neg_edged_vertices,
                            unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    bool no_vertex_removed = true;

    // vector<pair<int, Weight>> vec(pos_edged_vertices.begin(), pos_edged_vertices.end());
    vector<pair<int, double>> vec(pos_edged_vertices.begin(), pos_edged_vertices.end());
    sort(vec.begin(), vec.end(), compareByValue);
    for (auto u_it = vec.begin(); u_it != vec.end();)
    {
        int u = u_it->first;
        if (Check_Edge_vertex(u, graph, cc, k,
                              seq2att, Removed_vertices,
                              pos_edged_vertices, neg_edged_vertices, sim_graph))
        {
            no_vertex_removed = false;
            pos_edged_vertices.erase(u);
            // ns = (ns * cc.size()) / (cc.size() - 1);
            ns = (ns * (cc.size() + 1)) / cc.size();
            // cc.erase(find(cc.begin(), cc.end(), u));
            // ns = NSize(cc, seq2att);
            break;
        }
        else
        {
            u_it++;
        }
    }

    // for (auto u_it = pos_edged_vertices.begin(); u_it != pos_edged_vertices.end();)
    // {
    //     int u = *u_it;
    //     if (Check_Edge_vertex(u, graph, k,
    //                           seq2att, Removed_vertices,
    //                           pos_edged_vertices, neg_edged_vertices))
    //     {
    //         no_vertex_removed = false;
    //         u_it = pos_edged_vertices.erase(u_it);
    //         ns = (ns * cc.size()) / (cc.size() - 1);
    //         cc.erase(find(cc.begin(), cc.end(), u));
    //         // ns = NSize(cc, seq2att);
    //         break;
    //     }
    //     else
    //     {
    //         u_it++;
    //     }
    // }

    // else
    // {
    //     for (auto u_it = neg_edged_vertices.begin(); u_it != neg_edged_vertices.end();)
    //     {
    //         int u = *u_it;
    //         if (Check_Edge_vertex(u, graph, k,
    //                               seq2att, Removed_vertices,
    //                               pos_edged_vertices, neg_edged_vertices))
    //         {
    //             no_vertex_removed = false;
    //             u_it = neg_edged_vertices.erase(u_it);
    //             ns = (ns * cc.size() - 1) / (cc.size() - 1);
    //             cc.erase(find(cc.begin(), cc.end(), u));
    //             // ns = NSize(cc, seq2att);
    //             break;
    //         }
    //         else{
    //             u_it++;
    //         }
    //     }
    // }
    return no_vertex_removed;
}
bool First_Neg_Edged_Vertex(unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, int k, Weight &ns,
                            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, vector<VertexInfo> &Removed_vertices,
                            unordered_map<int, double> &pos_edged_vertices, unordered_map<int, double> &neg_edged_vertices,
                            unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    bool no_vertex_removed = true;
    // vector<pair<int, Weight>> vec(neg_edged_vertices.begin(), neg_edged_vertices.end());
    vector<pair<int, double>> vec(neg_edged_vertices.begin(), neg_edged_vertices.end());
    sort(vec.begin(), vec.end(), compareByValue);
    for (auto u_it = vec.begin(); u_it != vec.end();)
    {
        int u = u_it->first;
        if (Check_Edge_vertex(u, graph, cc, k,
                              seq2att, Removed_vertices,
                              pos_edged_vertices, neg_edged_vertices, sim_graph))
        {
            no_vertex_removed = false;
            neg_edged_vertices.erase(u);
            ns = (ns * (cc.size() + 1) - 1) / cc.size();
            // ns = (ns * cc.size() - 1) / (cc.size() - 1);
            // cc.erase(find(cc.begin(), cc.end(), u));
            // ns = NSize(cc, seq2att);
            break;
        }
        else
        {
            u_it++;
        }
    }

    // else
    // {
    //     for (auto u_it = pos_edged_vertices.begin(); u_it != pos_edged_vertices.end();)
    //     {
    //         int u = *u_it;
    //         if (Check_Edge_vertex(u, graph, k,
    //                               seq2att, Removed_vertices,
    //                               pos_edged_vertices, neg_edged_vertices))
    //         {
    //             no_vertex_removed = false;
    //             u_it = pos_edged_vertices.erase(u_it);
    //             ns = (ns * cc.size()) / (cc.size() - 1);
    //             cc.erase(find(cc.begin(), cc.end(), u));
    //             // ns = NSize(cc, seq2att);
    //             break;
    //         }
    //         else{
    //             u_it++;
    //         }
    //     }
    // }
    return no_vertex_removed;
}
// void printGraphAsEdgeList(const std::unordered_map<int, std::unordered_map<int, int>> &graph)
// {
//     std::cout << "edges = [\n";
//     for (const auto &[src, edges] : graph)
//     {
//         for (const auto &[dst, weight] : edges)
//         {
//             std::cout << "    (" << src << ", " << dst << ", " << weight << "),\n";
//         }
//     }
//     std::cout << "]\n";
// }
void printGraphAsEdgeList(const std::unordered_map<int, std::unordered_map<int, int>> &graph, vector<int> &cc)
{
    std::cout << "edges = [\n";
    for (int src : cc)
    {
        for (const auto &[dst, weight] : graph.at(src))
        {
            std::cout << "    (" << src << ", " << dst << ", " << weight << "),\n";
        }
    }
    std::cout << "]\n";
}
// void shrink(std::unordered_map<int, std::unordered_map<int, int>> &graph, vector<int> &cc, Weight &bottom_score,
//             map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
//             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//             unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     vector<VertexInfo> Removed_vertices;
//     // unordered_set<int> pos_edged_vertices;
//     // unordered_set<int> neg_edged_vertices;
//     unordered_map<int, Weight> pos_edged_vertices;
//     unordered_map<int, Weight> neg_edged_vertices;
//     for (int v : cc)
//     {
//         // if(!graph.count(src)) continue;

//         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//     }
//     // double ns = NSize(edge_truss, seq2att);
//     // double ns = 0.0;
//     // for (int u : cc)
//     // {
//     //     if (!seq2att[u].first)
//     //     {
//     //         ns++;
//     //     }
//     // }
//     // ns /= cc.size();
//     Weight ns = NSize(cc, seq2att);

//     while (cc.size() > k && (ns < thd1 || ns > thd2))
//     {
//         if (ns < thd1)
//         {
//             // remove the pos vertices
//             if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//                                        seq2att, Removed_vertices,
//                                        pos_edged_vertices, neg_edged_vertices, sim_graph))
//             {
//                 // break;
//                 if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//                                            seq2att, Removed_vertices,
//                                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//                 {
//                     break;
//                 }
//             }
//         }
//         else if (ns > thd2)
//         {
//             // remove the neg vertices
//             if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//                                        seq2att, Removed_vertices,
//                                        pos_edged_vertices, neg_edged_vertices, sim_graph))
//             {
//                 // break;
//                 if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//                                            seq2att, Removed_vertices,
//                                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//                 {
//                     break;
//                 }
//             }
//         }
//         // else
//         // {
//         //     Weight score = WEIGHT_MAX;
//         //     for (int i = 0; i < cc.size(); i++)
//         //     {
//         //         int u = cc[i];

//         //         for (int j = i + 1; j < cc.size(); j++)
//         //         {
//         //             int v = cc[j];

//         //             score = min(score, sim_graph[u][v]);
//         //         }
//         //     }
//         //     cout << "score: " << score.numerator << " / " << score.denominator << endl;
//         //     // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
//         //     if (score > bottom_score)
//         //     {
//         //         bottom_score = score;
//         //         if (bottom_score == Weight(1, 526))
//         //         // if (bottom_score == Weight(1, 192))
//         //         // if (bottom_score == Weight(1, 6336))
//         //         {
//         //             // printGraphAsEdgeList(graph, cc);
//         //             int a = 1;
//         //         }
//         //         // remove vertex with score less than bottom_score and maintain the truss
//         //         unordered_map<int, int> removedV_degree;
//         //         for (int u : cc)
//         //         {
//         //             if (removedV_degree.find(u) != removedV_degree.end())
//         //                 continue;
//         //             for (auto &v_s : sim_graph[u])
//         //             {
//         //                 int v = v_s.first;
//         //                 if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
//         //                 {
//         //                     removedV_degree.emplace(u, graph[u].size());
//         //                     removedV_degree.emplace(v, graph[v].size());
//         //                 }
//         //             }
//         //         }
//         //         // 按度升序排列，优先删除度小的，在删除过程中更新要删除的点及其度
//         //         // vector<int> sorted_elbys(removedV_degree.size());
//         //         // map<int, int> sorted_ep;
//         //         // unordered_map<int, int> svp;
//         //         // vector<int> bucket;
//         //         // int degmax = 0;
//         //         // for(int u : cc){
//         //         //     int deg = graph[u].size();
//         //         //     degmax = max(degmax, deg);
//         //         // }
//         //         // bucket.resize(degmax + 1, 0);
//         //         // for (auto &u_d : removedV_degree)
//         //         // {
//         //         //     bucket[u_d.second]++; // bucket[i] = the number of edges whose support = i, 同一支持度的边有多少
//         //         // }
//         //         // int p = 0;
//         //         // for (int j = 0; j < degmax + 1; j++)
//         //         // {
//         //         //     int tmp = bucket[j];
//         //         //     bucket[j] = p; // record the first index of edge whose support=j after sort
//         //         //     p += tmp;
//         //         // }
//         //         // for (auto &u_d : removedV_degree)
//         //         // {
//         //         //     sorted_elbys[bucket[u_d.second]] = u_d.first;     // edges sort as support increasing
//         //         //     sorted_ep.emplace(u_d.first, bucket[u_d.second]); // edge : index after sort
//         //         //     svp.emplace(u_d.second, bucket[u_d.second]);      // support: the first index of edge whose support=support after sort
//         //         //     bucket[u_d.second]++;
//         //         // }
//         //         // vector<VertexInfo> changed_edges;
//         //         // set<int> del_vertex;
//         //         // for (int i = 0; i < sorted_elbys.size(); i++){
//         //         //     int u = sorted_elbys[i];
//         //         //     if(del_vertex.count(u)){
//         //         //         continue;
//         //         //     }
//         //         //     set<int> remove_vertex;
//         //         //     for (auto &v_s : sim_graph[u])
//         //         //     {
//         //         //         int v = v_s.first;
//         //         //         if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
//         //         //         {
//         //         //             changed_edges.push_back(MaintainTruss_DeleteVertex(u, graph, k, remove_vertex));
//         //         //             break;
//         //         //         }
//         //         //     }
//         //         // }

//         //         vector<VertexInfo> changed_edges;
//         //         vector<set<int>> del_vertices;

//         //         while (!removedV_degree.empty())
//         //         {
//         //             // find the vertex with the smallest degree
//         //             int minV = -1;
//         //             int minD = INT_MAX;
//         //             for (auto &it : removedV_degree)
//         //             {
//         //                 int v = it.first;
//         //                 int d = it.second;
//         //                 if (d < minD)
//         //                 {
//         //                     minV = v;
//         //                     minD = d;
//         //                 }
//         //             }
//         //             set<int> del_vertex;
//         //             changed_edges.push_back(MaintainTruss_DeleteVertex(minV, graph, k, del_vertex));
//         //             del_vertices.push_back(del_vertex);
//         //             // update cc
//         //             cc.erase(find(cc.begin(), cc.end(), minV));
//         //             // update remove_degree
//         //             removedV_degree.erase(minV);
//         //             for (int u : del_vertex)
//         //             {
//         //                 cc.erase(find(cc.begin(), cc.end(), u));
//         //                 removedV_degree.erase(u);
//         //             }
//         //             for (auto it = removedV_degree.begin(); it != removedV_degree.end(); ) {
//         //                 int u = it->first;
//         //                 bool flag = false;
//         //                 for (auto &v_s : sim_graph[u])
//         //                 {
//         //                     int v = v_s.first;
//         //                     if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
//         //                     {
//         //                         flag = true;
//         //                         break;
//         //                     }
//         //                 }
//         //                 if (!flag)
//         //                 {
//         //                     it = removedV_degree.erase(it);
//         //                 }
//         //                 else
//         //                 {
//         //                     removedV_degree[u] = graph[u].size();
//         //                     it++;
//         //                 }
//         //             }
//         //             // // 判断是否还需要删除
//         //             // if (!del_vertex.count(minV) && cover_less_edge())
//         //             // {
//         //             //     changed_edges.push_back(MaintainTruss_DeleteVertex(u, graph, k, del_vertex));
//         //             // }
//         //             // vec.pop_back();
//         //         }
//         //         // find connected components
//         //         vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
//         //         // calculate all added vertex-pair score
//         //         // set<int> added_vertices;
//         //         for (auto &cc : connectedComponents)
//         //         {
//         //             if (cc.size() <= 1)
//         //                 continue;
//         //             // repeat the process
//         //             shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//         //         }

//         //         // for (int u : removed_vertices)
//         //         // {
//         //         //     if (del_vertex.count(u))
//         //         //         continue;
//         //         //     changed_edges.push_back(MaintainTruss_DeleteVertex(u, graph, k, del_vertex));
//         //         //     // del_vertex.insert(u);
//         //         // }
//         //         for (auto &it : del_vertices)
//         //         {
//         //             for (int v : it)
//         //             {
//         //                 graph.emplace(v, unordered_map<int, int>());
//         //             }
//         //         }

//         //         while (!changed_edges.empty())
//         //         {
//         //             auto &it = changed_edges.back();

//         //             recoverTruss(it, graph);
//         //             changed_edges.pop_back();
//         //         }
//         //     }

//         //     if (First_Neg_Edged_Vertex(graph, k, ns, cc,
//         //                                seq2att, Removed_vertices,
//         //                                pos_edged_vertices, neg_edged_vertices, sim_graph))
//         //     {
//         //         if (First_Pos_Edged_Vertex(graph, k, ns, cc,
//         //                                    seq2att, Removed_vertices,
//         //                                    pos_edged_vertices, neg_edged_vertices, sim_graph))
//         //         {
//         //             break;
//         //         }
//         //     }
//         // }
//     }
//     if (ns >= thd1 && ns <= thd2)
//     {
//         Weight score = WEIGHT_MAX;
//         for (int i = 0; i < cc.size(); i++)
//         {
//             int u = cc[i];

//             for (int j = i + 1; j < cc.size(); j++)
//             {
//                 int v = cc[j];

//                 score = min(score, sim_graph[u][v]);
//             }
//         }
//         // cout << "score: " << score.numerator << " / " << score.denominator << endl;
//         // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
//         if (score > bottom_score)
//         {
//             bottom_score = score;
//             if (bottom_score == Weight(1, 526))
//             // if (bottom_score == Weight(1, 192))
//             // if (bottom_score == Weight(1, 6336))
//             {
//                 // printGraphAsEdgeList(graph, cc);
//                 int a = 1;
//             }
//             // printGraphAsEdgeList(graph, cc);
//             // remove vertex with score less than bottom_score and maintain the truss
//             unordered_map<int, int> removedV_degree;
//             for (int u : cc)
//             {
//                 if (removedV_degree.find(u) != removedV_degree.end())
//                     continue;
//                 for (auto &v_s : sim_graph[u])
//                 {
//                     int v = v_s.first;
//                     if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
//                     {
//                         removedV_degree.emplace(u, graph[u].size());
//                         removedV_degree.emplace(v, graph[v].size());
//                     }
//                 }
//             }

//             vector<VertexInfo> changed_edges;
//             vector<set<int>> del_vertices;

//             while (!removedV_degree.empty())
//             {
//                 // find the vertex with the smallest degree
//                 int minV = -1;
//                 int minD = INT_MAX;
//                 for (auto &it : removedV_degree)
//                 {
//                     int v = it.first;
//                     int d = it.second;
//                     if (d < minD)
//                     {
//                         minV = v;
//                         minD = d;
//                     }
//                 }
//                 set<int> del_vertex;
//                 changed_edges.push_back(MaintainTruss_DeleteVertex(minV, graph, k, del_vertex));
//                 del_vertices.push_back(del_vertex);
//                 // update cc
//                 cc.erase(find(cc.begin(), cc.end(), minV));
//                 // update remove_degree
//                 removedV_degree.erase(minV);
//                 for (int u : del_vertex)
//                 {
//                     cc.erase(find(cc.begin(), cc.end(), u));
//                     removedV_degree.erase(u);
//                 }
//                 for (auto it = removedV_degree.begin(); it != removedV_degree.end();)
//                 {
//                     int u = it->first;
//                     bool flag = false;
//                     for (auto &v_s : sim_graph[u])
//                     {
//                         int v = v_s.first;
//                         if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
//                         {
//                             flag = true;
//                             break;
//                         }
//                     }
//                     if (!flag)
//                     {
//                         it = removedV_degree.erase(it);
//                     }
//                     else
//                     {
//                         removedV_degree[u] = graph[u].size();
//                         it++;
//                     }
//                 }
//             }

//             // find connected components
//             vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
//             // calculate all added vertex-pair score
//             // set<int> added_vertices;
//             for (auto &cc : connectedComponents)
//             {
//                 if (cc.size() <= 1)
//                     continue;
//                 // repeat the process
//                 // printGraphAsEdgeList(graph, cc);
//                 shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//             }
//             for (auto &it : del_vertices)
//             {
//                 for (int v : it)
//                 {
//                     graph.emplace(v, unordered_map<int, int>());
//                 }
//             }

//             while (!changed_edges.empty())
//             {
//                 auto &it = changed_edges.back();

//                 recoverTruss(it, graph);
//                 changed_edges.pop_back();
//             }
//         }
//     }

//     // recovery
//     // if (score_id == score_edges.size())
//     //     break;
//     while (!Removed_vertices.empty())
//     {
//         VertexInfo v = Removed_vertices.back();
//         RecoverEdges(v.vertex, v.neighbors, v.changed_edges, graph);
//         Removed_vertices.pop_back();
//     }
// }
// void shrink(std::unordered_map<int, std::unordered_map<int, int>> &graph, vector<int> &cc, Weight &bottom_score,
//             map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
//             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//             unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     vector<VertexInfo> Removed_vertices;
//     // unordered_set<int> pos_edged_vertices;
//     // unordered_set<int> neg_edged_vertices;
//     unordered_map<int, double> pos_edged_vertices;
//     unordered_map<int, double> neg_edged_vertices;
//     for (int v : cc)
//     {
//         // if(!graph.count(src)) continue;

//         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//     }
//     // double ns = NSize(edge_truss, seq2att);
//     // double ns = 0.0;
//     // for (int u : cc)
//     // {
//     //     if (!seq2att[u].first)
//     //     {
//     //         ns++;
//     //     }
//     // }
//     // ns /= cc.size();
//     Weight ns = NSize(cc, seq2att);

//     while (cc.size() > k && (ns < thd1 || ns > thd2))
//     {
//         bool has_vertex = true;
//         if (ns < thd1)
//         {
//             // remove the pos vertices
//             if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//                                        seq2att, Removed_vertices,
//                                        pos_edged_vertices, neg_edged_vertices, sim_graph))
//             {
//                 break;
//                 // if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//                 //                            seq2att, Removed_vertices,
//                 //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//                 // {
//                 //     // has_vertex =  false;
//                 //     break;
//                 // }
//             }
//         }
//         else if (ns > thd2)
//         {
//             // remove the neg vertices
//             if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//                                        seq2att, Removed_vertices,
//                                        pos_edged_vertices, neg_edged_vertices, sim_graph))
//             {
//                 break;
//                 // if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//                 //                            seq2att, Removed_vertices,
//                 //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//                 // {
//                 //     // has_vertex =  false;
//                 //     break;
//                 // }
//             }
//         }
//     }
//     if (ns >= thd1 && ns <= thd2)
//     {
//         Weight score = WEIGHT_MAX;
//         for (int i = 0; i < cc.size(); i++)
//         {
//             int u = cc[i];

//             for (int j = i + 1; j < cc.size(); j++)
//             {
//                 int v = cc[j];

//                 score = min(score, sim_graph[u][v]);
//             }
//         }
//         // cout << "score: " << score.numerator << " / " << score.denominator << endl;
//         // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
//         if (score > bottom_score)
//         {
//             bottom_score = score;
//             if (bottom_score == Weight(1, 526))
//             // if (bottom_score == Weight(1, 192))
//             // if (bottom_score == Weight(1, 6336))
//             {
//                 // printGraphAsEdgeList(graph, cc);
//                 int a = 1;
//             }
//             printGraphAsEdgeList(graph, cc);
//             // remove vertex with score less than bottom_score and maintain the truss
//             unordered_map<int, int> removedV_degree;
//             for (int u : cc)
//             {
//                 if (removedV_degree.find(u) != removedV_degree.end())
//                     continue;
//                 for (auto &v_s : sim_graph[u])
//                 {
//                     int v = v_s.first;
//                     if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
//                     {
//                         removedV_degree.emplace(u, graph[u].size());
//                         removedV_degree.emplace(v, graph[v].size());
//                     }
//                 }
//             }

//             vector<VertexInfo> changed_edges;
//             vector<set<int>> del_vertices;

//             while (!removedV_degree.empty())
//             {
//                 // find the vertex with the smallest degree
//                 int minV = -1;
//                 int minD = INT_MAX;
//                 for (auto &it : removedV_degree)
//                 {
//                     int v = it.first;
//                     int d = it.second;
//                     if (d < minD)
//                     {
//                         minV = v;
//                         minD = d;
//                     }
//                 }
//                 set<int> del_vertex;
//                 changed_edges.push_back(MaintainTruss_DeleteVertex(minV, graph, k, del_vertex));
//                 del_vertices.push_back(del_vertex);
//                 // update cc
//                 cc.erase(find(cc.begin(), cc.end(), minV));
//                 // update remove_degree
//                 removedV_degree.erase(minV);
//                 for (int u : del_vertex)
//                 {
//                     cc.erase(find(cc.begin(), cc.end(), u));
//                     removedV_degree.erase(u);
//                 }
//                 for (auto it = removedV_degree.begin(); it != removedV_degree.end();)
//                 {
//                     int u = it->first;
//                     bool flag = false;
//                     for (auto &v_s : sim_graph[u])
//                     {
//                         int v = v_s.first;
//                         if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
//                         {
//                             flag = true;
//                             break;
//                         }
//                     }
//                     if (!flag)
//                     {
//                         it = removedV_degree.erase(it);
//                     }
//                     else
//                     {
//                         removedV_degree[u] = graph[u].size();
//                         it++;
//                     }
//                 }
//             }

//             // find connected components
//             vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
//             // calculate all added vertex-pair score
//             // set<int> added_vertices;
//             for (auto &cc : connectedComponents)
//             {
//                 if (cc.size() <= 1)
//                     continue;
//                 // repeat the process
//                 // printGraphAsEdgeList(graph, cc);
//                 cout << "delete edges" << endl;
//                 shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//             }
//             for (auto &it : del_vertices)
//             {
//                 for (int v : it)
//                 {
//                     graph.emplace(v, unordered_map<int, int>());
//                 }
//             }

//             while (!changed_edges.empty())
//             {
//                 auto &it = changed_edges.back();

//                 recoverTruss(it, graph);
//                 changed_edges.pop_back();
//             }
//         }
//     }

//     // recovery
//     // if (score_id == score_edges.size())
//     //     break;
//     while (!Removed_vertices.empty())
//     {
//         VertexInfo v = Removed_vertices.back();
//         RecoverEdges(v.vertex, v.neighbors, v.changed_edges, graph);
//         Removed_vertices.pop_back();
//     }
// }
void greedy_remove(vector<VertexInfo> &changed_edges,
                   vector<set<int>> &del_vertices, vector<int> &cc, unordered_map<int, unordered_map<int, int>> &graph,
                   unordered_map<int, unordered_map<int, Weight>> &sim_graph,
                   Weight bottom_score, int k)
{
    unordered_map<int, unordered_set<int>> removedV2lessNum;
    for (int u : cc)
    {
        if (removedV2lessNum.find(u) != removedV2lessNum.end())
            continue;

        // int num = 0;
        for (auto &v_s : sim_graph[u])
        {
            int v = v_s.first;
            if (u < v)
            {
                continue;
            }
            if (count(cc.begin(), cc.end(), v) && v_s.second <= bottom_score)
            {
                if (!removedV2lessNum.count(u))
                {
                    removedV2lessNum.emplace(u, unordered_set<int>());
                }
                removedV2lessNum[u].insert(v);
                if (!removedV2lessNum.count(v))
                {
                    removedV2lessNum.emplace(v, unordered_set<int>());
                }
                removedV2lessNum[v].insert(u);
            }
        }
    }

    // vector<VertexInfo> changed_edges;
    // vector<set<int>> del_vertices;

    while (!removedV2lessNum.empty())
    {
        // find the vertex with the smallest degree
        int maxV = -1;
        int maxN = 0;
        for (auto &it : removedV2lessNum)
        {
            int v = it.first;
            int d = it.second.size();
            if (d > maxN)
            {
                maxV = v;
                maxN = d;
            }
        }
        set<int> del_vertex;

        changed_edges.push_back(MaintainTruss_DeleteVertex(maxV, graph, k, del_vertex));
        del_vertices.push_back(del_vertex);
        // update cc
        cc.erase(find(cc.begin(), cc.end(), maxV));
        // update remove_degree
        for (auto &vs : removedV2lessNum[maxV])
        {
            removedV2lessNum[vs].erase(maxV);
            if (removedV2lessNum[vs].empty())
            {
                removedV2lessNum.erase(vs);
            }
        }
        removedV2lessNum.erase(maxV);
        for (int u : del_vertex)
        {
            cc.erase(find(cc.begin(), cc.end(), u));
            for (auto &vs : removedV2lessNum[u])
            {
                removedV2lessNum[vs].erase(u);
                if (removedV2lessNum[vs].empty())
                {
                    removedV2lessNum.erase(vs);
                }
            }

            removedV2lessNum.erase(u);
        }
    }
}
// void shrink(std::unordered_map<int, std::unordered_map<int, int>> &graph, vector<int> &cc, Weight &bottom_score,
//             map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
//             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//             unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     vector<VertexInfo> Removed_vertices;
//     // unordered_set<int> pos_edged_vertices;
//     // unordered_set<int> neg_edged_vertices;
//     unordered_map<int, double> pos_edged_vertices;
//     unordered_map<int, double> neg_edged_vertices;
//     for (int v : cc)
//     {
//         // if(!graph.count(src)) continue;

//         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//     }
//     // double ns = NSize(edge_truss, seq2att);
//     // double ns = 0.0;
//     // for (int u : cc)
//     // {
//     //     if (!seq2att[u].first)
//     //     {
//     //         ns++;
//     //     }
//     // }
//     // ns /= cc.size();
//     Weight ns = NSize(cc, seq2att);

//     while (cc.size() > k && (ns < thd1 || ns > thd2))
//     {
//         bool has_vertex = true;
//         if (ns < thd1)
//         {
//             // remove the pos vertices
//             if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//                                        seq2att, Removed_vertices,
//                                        pos_edged_vertices, neg_edged_vertices, sim_graph))
//             {
//                 break;
//                 // if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//                 //                            seq2att, Removed_vertices,
//                 //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//                 // {
//                 //     // has_vertex =  false;
//                 //     break;
//                 // }
//             }
//         }
//         else if (ns > thd2)
//         {
//             // remove the neg vertices
//             if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//                                        seq2att, Removed_vertices,
//                                        pos_edged_vertices, neg_edged_vertices, sim_graph))
//             {
//                 printGraphAsEdgeList(graph,cc);
//                 break;
//                 // if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//                 //                            seq2att, Removed_vertices,
//                 //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//                 // {
//                 //     // has_vertex =  false;
//                 //     break;
//                 // }
//             }
//         }
//     }
//     if (ns >= thd1 && ns <= thd2)
//     {
//         Weight score = WEIGHT_MAX;
//         for (int i = 0; i < cc.size(); i++)
//         {
//             int u = cc[i];

//             for (int j = i + 1; j < cc.size(); j++)
//             {
//                 int v = cc[j];

//                 score = min(score, sim_graph[u][v]);
//             }
//         }
//         // cout << "score: " << score.numerator << " / " << score.denominator << endl;
//         // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
//         if (score > bottom_score)
//         {
//             bottom_score = score;
//             if (bottom_score == Weight(1, 526))
//             // if (bottom_score == Weight(1, 192))
//             // if (bottom_score == Weight(1, 6336))
//             {
//                 // printGraphAsEdgeList(graph, cc);
//                 int a = 1;
//             }
//             printGraphAsEdgeList(graph, cc);
//             // remove vertex with score less than bottom_score and maintain the truss
//             vector<VertexInfo> changed_edges;
//             vector<set<int>> del_vertices;
//             greedy_remove(changed_edges, del_vertices, cc, graph, sim_graph, bottom_score, k);
//             // find connected components
//             vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
//             // calculate all added vertex-pair score
//             // set<int> added_vertices;
//             for (auto &cc : connectedComponents)
//             {
//                 if (cc.size() <= 1)
//                     continue;
//                 // repeat the process
//                 // printGraphAsEdgeList(graph, cc);
//                 cout << "delete edges" << endl;
//                 shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//             }
//             for (auto &it : del_vertices)
//             {
//                 for (int v : it)
//                 {
//                     graph.emplace(v, unordered_map<int, int>());
//                 }
//             }

//             while (!changed_edges.empty())
//             {
//                 auto &it = changed_edges.back();

//                 recoverTruss(it, graph);
//                 changed_edges.pop_back();
//             }
//         }
//     }

//     // recovery
//     // if (score_id == score_edges.size())
//     //     break;
//     while (!Removed_vertices.empty())
//     {
//         VertexInfo v = Removed_vertices.back();
//         RecoverEdges(v.vertex, v.neighbors, v.changed_edges, graph);
//         Removed_vertices.pop_back();
//     }
// }
void BCCUtil(unordered_map<int, unordered_map<int, int>> &graph, set<int> &cv, int u, unordered_map<int, int> &disc, unordered_map<int, int> &low, vector<int> &st,
             unordered_map<int, int> &parent, int &time)
{
    disc[u] = low[u] = ++time;
    int children = 0;
    st.push_back(u);
    for (auto &v_sup : graph[u])
    {
        int v = v_sup.first;
        if (!disc.count(v))
            continue;
        if (disc[v] == -1)
        {
            children++;
            parent[v] = u;

            BCCUtil(graph, cv, v, disc, low, st, parent, time);

            low[u] = min(low[u], low[v]);

            if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && low[v] >= disc[u]))
            {
                cv.insert(u);
                vector<int> bcc;
                int x;
                do
                {
                    x = st.back();
                    bcc.push_back(x);
                    st.pop_back();
                } while (x != v);
            }
        }

        else if (v != parent[u])
        {
            low[u] = min(low[u], disc[v]);
        }
    }
}
void BCC(unordered_map<int, unordered_map<int, int>> &graph, set<int> &cv)
{
    unordered_map<int, int> disc;
    unordered_map<int, int> low;
    unordered_map<int, int> parent;
    int time = 0;
    int index = 0;
    // vector<vector<int>> bccs;
    vector<int> st;

    // get neighbors of delete vertex
    //  int dv;
    // vector<int> trueSub;
    // Initialize disc and low, and parent arrays
    for (auto &u_v_sup : graph)
    {
        int i = u_v_sup.first;
        disc[i] = NIL;
        low[i] = NIL;
        parent[i] = NIL;
        // trueSub.push_back(i);
    }

    for (auto &u_v_sup : graph)
    {
        int i = u_v_sup.first;
        auto t = disc.find(i);
        if (t != disc.end() && t->second == NIL)
            BCCUtil(graph, cv, i, disc, low, st, parent, time);

        // If stack is not empty, pop all edges from stack
        while (st.size() > 0)
        {

            // bcc.push_back(st.back());
            st.pop_back();
        }
    }
}
void dijkstra(int sourceNode, int targetNode,
              unordered_map<int, unordered_map<int, double>> &graph, vector<int> &path)
{
    // cout << "dijkstra" << endl;
    unordered_map<int, int> predecessor;
    // 使用unordered_map动态记录距离
    unordered_map<int, double> distance;
    unordered_map<int, bool> visited;

    for (auto &e : graph)
    {
        int node = e.first;
        distance[node] = numeric_limits<double>::max();
        visited[node] = false;
    }

    priority_queue<Node, vector<Node>, greater<Node>> pq;

    // 初始化源点
    distance[sourceNode] = 0;
    pq.push({sourceNode, 0});

    while (!pq.empty())
    {
        Node current = pq.top();
        pq.pop();

        int u = current.vertex;

        // 如果当前节点已访问过，跳过
        if (visited[u])
            continue;

        visited[u] = true;

        // 如果到达目标节点，返回最短距离
        if (u == targetNode)
        {
            // return distance[u];
            break;
        }

        // 遍历当前节点的邻居
        if (graph.find(u) != graph.end())
        {
            for (const auto &edge : graph.at(u))
            {
                int v = edge.first;
                double weight = edge.second;

                // 如果找到更短的路径，则更新距离和前驱节点，并将节点加入队列
                if (!visited[v] && distance[u] + weight < distance[v])
                {
                    distance[v] = distance[u] + weight;
                    predecessor[v] = u; // 记录前驱
                    pq.push({v, distance[v]});
                }
            }
        }
    }

    // 如果目标节点不可达，返回 -1 表示无解
    // return -1;
    int current = targetNode;

    // 回溯前驱数组，构建路径
    while (predecessor.find(current) != predecessor.end())
    {
        if (current != targetNode)
            path.push_back(current);
        current = predecessor.at(current);
    }
}
// void shrink(std::unordered_map<int, std::unordered_map<int, int>> &graph, vector<int> &cc, Weight &bottom_score,
//             map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
//             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//             unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     Weight ns = NSize(cc, seq2att);
//     if (cc.size() == k)
//     {
//         if (ns >= thd1 && ns <= thd2)
//         {
//             Weight score = WEIGHT_MAX;
//             for (int i = 0; i < cc.size(); i++)
//             {
//                 int u = cc[i];

//                 for (int j = i + 1; j < cc.size(); j++)
//                 {
//                     int v = cc[j];

//                     score = min(score, sim_graph[u][v]);
//                 }
//             }
//             bottom_score = score;
//             // printGraphAsEdgeList(graph, cc);
//         }
//         return;
//     }
//     if (thd2 < 1 && ns == Weight(1, 1))
//     {
//         return;
//     }
//     if (thd1 > 0 && ns == Weight(0, 1))
//     {
//         return;
//     }
//     vector<VertexInfo> Removed_vertices;
//     // unordered_set<int> pos_edged_vertices;
//     // unordered_set<int> neg_edged_vertices;
//     unordered_map<int, double> pos_edged_vertices;
//     unordered_map<int, double> neg_edged_vertices;
//     for (int v : cc)
//     {
//         // if(!graph.count(src)) continue;
//         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//     }
//     // double ns = NSize(edge_truss, seq2att);
//     // double ns = 0.0;
//     // for (int u : cc)
//     // {
//     //     if (!seq2att[u].first)
//     //     {
//     //         ns++;
//     //     }
//     // }
//     // ns /= cc.size();
//     while (cc.size() > k && (ns < thd1 || ns > thd2))
//     {
//         // bool has_vertex = true;

//         if (ns < thd1)
//         {
//             // remove the pos vertices
//             vector<pair<int, double>> vec(pos_edged_vertices.begin(), pos_edged_vertices.end());
//             sort(vec.begin(), vec.end(), compareByValue);
//             auto u_it = vec.begin();
//             for (; u_it != vec.end();)
//             {
//                 int u = u_it->first;
//                 unordered_set<int> deleted_vertices;
//                 unordered_map<int, unordered_map<int, int>> deleted_edges;
//                 // vector<Edge> changed_edges;
//                 Deleted_vertex(u, graph, cc, k, deleted_vertices, deleted_edges);
//                 if (deleted_vertices.empty())
//                 {
//                     u_it++;
//                     continue;
//                 }
//                 int neg_num = 0;
//                 for (int v : deleted_vertices)
//                 {
//                     if (!seq2att[v].first)
//                     {
//                         neg_num++;
//                     }
//                 }
//                 Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - deleted_vertices.size());
//                 if (new_ns < ns)
//                 {
//                     u_it++;
//                     pos_edged_vertices.erase(u);
//                     continue;
//                 }
//                 else
//                 {
//                     ns = new_ns;
//                     // update pos/neg edged vertices
//                     for (int u : deleted_vertices)
//                     {
//                         if (pos_edged_vertices.count(u))
//                         {
//                             pos_edged_vertices.erase(u);
//                         }
//                         else if (neg_edged_vertices.count(u))
//                         {
//                             neg_edged_vertices.erase(u);
//                         }
//                     }
//                     // unordered_map<int, double> pos_neighbors2delete_weight;
//                     // unordered_map<int, int> pos_neighbors2delete_count;
//                     // unordered_map<int, double> neg_neighbors2delete_weight;
//                     // unordered_map<int, int> neg_neighbors2delete_count;
//                     // for (int v : deleted_vertices)
//                     // {
//                     //     for (auto &w_s : graph[v])
//                     //     {
//                     //         int w = w_s.first;
//                     //         if (pos_edged_vertices.count(w))
//                     //         {
//                     //             Weight a = sim_graph[v][w] + 1;
//                     //             // int neighbor_count = graph[w].size();
//                     //             // double b = pos_edged_vertices[w] * neighbor_count;
//                     //             if(!pos_neighbors2delete_weight.count(w)){
//                     //                 pos_neighbors2delete_weight.emplace(w, 0.0);
//                     //                 pos_neighbors2delete_count.emplace(w, 0);
//                     //             }
//                     //             pos_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                     //             pos_neighbors2delete_count[w] += 1;
//                     //         }
//                     //         else if(neg_edged_vertices.count(w)){
//                     //             Weight a = sim_graph[v][w] + 1;
//                     //             // neg_edged_vertices[w] -= log(a.numerator) - log(a.denominator);
//                     //             if(!neg_neighbors2delete_weight.count(w)){
//                     //                 neg_neighbors2delete_weight.emplace(w, 0.0);
//                     //                 neg_neighbors2delete_count.emplace(w, 0);
//                     //             }
//                     //             neg_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                     //             neg_neighbors2delete_count[w] += 1;
//                     //         }
//                     //     }
//                     // }
//                     // for(auto &v_w : pos_neighbors2delete_weight){
//                     //     int v = v_w.first;
//                     //     double ori_weight = pos_edged_vertices[v] * graph[v].size();
//                     //     pos_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - pos_neighbors2delete_count[v]);
//                     // }
//                     // for(auto &v_w : neg_neighbors2delete_weight){
//                     //     int v = v_w.first;
//                     //     double ori_weight = neg_edged_vertices[v] * graph[v].size();
//                     //     neg_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - neg_neighbors2delete_count[v]);

//                     // }
//                     for (auto &v_s : pos_edged_vertices)
//                     {
//                         int v = v_s.first;
//                         double ori_w = pos_edged_vertices[v] * cc.size();
//                         for (int w : deleted_vertices)
//                         {
//                             Weight a = sim_graph[w][v] + 1;
//                             ori_w -= log(a.numerator) - log(a.denominator);
//                         }
//                         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                     }
//                     for (auto &v_s : neg_edged_vertices)
//                     {
//                         int v = v_s.first;
//                         double ori_w = neg_edged_vertices[v] * cc.size();
//                         for (int w : deleted_vertices)
//                         {
//                             Weight a = sim_graph[w][v] + 1;
//                             ori_w -= log(a.numerator) - log(a.denominator);
//                         }
//                         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                     }
//                     // for (auto &v_s : pos_edged_vertices)
//                     // {
//                     //     int v = v_s.first;
//                     //     for (int w : deleted_vertices)
//                     //     {
//                     //         Weight a = sim_graph[w][v] + 1;
//                     //         v_s.second -= log(a.numerator) - log(a.denominator);
//                     //     }
//                     // }
//                     // for (auto &v_s : neg_edged_vertices)
//                     // {
//                     //     int v = v_s.first;
//                     //     for (int w : deleted_vertices)
//                     //     {
//                     //         Weight a = sim_graph[w][v] + 1;
//                     //         v_s.second -= log(a.numerator) - log(a.denominator);
//                     //     }
//                     // }

//                     unordered_set<int> new_vertex;
//                     for (int u : deleted_vertices)
//                     {
//                         for (auto &v_sup : graph[u])
//                         {
//                             int v = v_sup.first;
//                             // graph[v].erase(u); // remove u from v's neighbors

//                             // check whether has new edged vertex
//                             if (!new_vertex.count(v) && !pos_edged_vertices.count(v) && !neg_edged_vertices.count(v))
//                             {
//                                 new_vertex.insert(v);
//                                 // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                             }
//                         }
//                     }
//                     // update graph
//                     for (int u : deleted_vertices)
//                     {
//                         unordered_map<int, int> neighbors = graph[u];
//                         vector<Edge> changed_edges;
//                         for (auto &v_sup : neighbors)
//                         {
//                             int v = v_sup.first;
//                             graph[v].erase(u);
//                             graph[u].erase(v);
//                             for (auto &it1 : graph[v])
//                             {
//                                 int w = it1.first;
//                                 if (graph[u].count(w))
//                                 {
//                                     --graph[v][w];
//                                     --graph[w][v];
//                                     changed_edges.push_back(pair<int, int>(v, w));
//                                 }
//                             }
//                         }

//                         cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
//                         graph.erase(u);
//                         VertexInfo removeV;
//                         removeV.vertex = u;
//                         removeV.neighbors = move(neighbors);
//                         removeV.changed_edges = move(changed_edges);
//                         Removed_vertices.push_back(removeV);
//                     }
//                     // printGraphAsEdgeList(graph, cc);
//                     // for (auto &u_v : deleted_edges)
//                     // {
//                     //     int u = u_v.first;
//                     //     for (auto &v_sup : u_v.second)
//                     //     {
//                     //         int v = v_sup.first;
//                     //         graph[u].erase(v);
//                     //     }
//                     // }

//                     for (auto &v : new_vertex)
//                     {
//                         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                     }

//                     // for (auto &u_v : changed_edges)
//                     // {
//                     //     int u = u_v.first;
//                     //     int v = u_v.second;
//                     //     graph[u][v]--;
//                     //     graph[v][u]--;
//                     // }

//                     break;
//                 }
//             }
//             if (u_it == vec.end())
//             {
//                 break;
//             }
//         }
//         else if (ns > thd2)
//         {
//             unordered_map<int, double> new_neg_edged_vertices;
//             for(auto &u_w : neg_edged_vertices){
//                 int u = u_w.first;
//                 int pos_neighbor_size = 0;
//                 for(auto &u_v : graph[u]){
//                     int v = u_v.first;
//                     if(seq2att[v].first){
//                         pos_neighbor_size++;
//                     }
//                 }
//                 new_neg_edged_vertices.emplace(u, u_w.second + pos_neighbor_size);
//             }
//             // vector<pair<int, double>> vec(neg_edged_vertices.begin(), neg_edged_vertices.end());
//             vector<pair<int, double>> vec(new_neg_edged_vertices.begin(), new_neg_edged_vertices.end());
//             sort(vec.begin(), vec.end(), compareByValue);
//             auto u_it = vec.begin();
//             for (; u_it != vec.end();)
//             {
//                 int u = u_it->first;
//                 if (u == 12)
//                 {
//                     int a = 1;
//                 }
//                 unordered_set<int> deleted_vertices;
//                 // unordered_map<int, unordered_map<int, int>> deleted_edges;
//                 unordered_map<int, unordered_map<int, int>> deleted_edges;
//                 Deleted_vertex(u, graph, cc, k, deleted_vertices, deleted_edges);
//                 if (deleted_vertices.empty())
//                 {
//                     u_it++;
//                     continue;
//                 }
//                 int neg_num = 0;
//                 for (int v : deleted_vertices)
//                 {
//                     if (!seq2att[v].first)
//                     {
//                         neg_num++;
//                     }
//                 }
//                 Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - deleted_vertices.size());
//                 if (new_ns > ns)
//                 {
//                     u_it++;
//                     neg_edged_vertices.erase(u);
//                     continue;
//                 }
//                 else
//                 {
//                     ns = new_ns;
//                     // update pos/neg edged vertices
//                     for (int u : deleted_vertices)
//                     {
//                         if (pos_edged_vertices.count(u))
//                         {
//                             pos_edged_vertices.erase(u);
//                         }
//                         else if (neg_edged_vertices.count(u))
//                         {
//                             neg_edged_vertices.erase(u);
//                         }
//                     }
//                     // unordered_map<int, double> pos_neighbors2delete_weight;
//                     // unordered_map<int, int> pos_neighbors2delete_count;
//                     // unordered_map<int, double> neg_neighbors2delete_weight;
//                     // unordered_map<int, int> neg_neighbors2delete_count;
//                     // for (int v : deleted_vertices)
//                     // {
//                     //     for (auto &w_s : graph[v])
//                     //     {
//                     //         int w = w_s.first;
//                     //         if (pos_edged_vertices.count(w))
//                     //         {
//                     //             Weight a = sim_graph[v][w] + 1;
//                     //             // int neighbor_count = graph[w].size();
//                     //             // double b = pos_edged_vertices[w] * neighbor_count;
//                     //             if(!pos_neighbors2delete_weight.count(w)){
//                     //                 pos_neighbors2delete_weight.emplace(w, 0.0);
//                     //                 pos_neighbors2delete_count.emplace(w, 0);
//                     //             }
//                     //             pos_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                     //             pos_neighbors2delete_count[w] += 1;
//                     //         }
//                     //         else if(neg_edged_vertices.count(w)){
//                     //             Weight a = sim_graph[v][w] + 1;
//                     //             // neg_edged_vertices[w] -= log(a.numerator) - log(a.denominator);
//                     //             if(!neg_neighbors2delete_weight.count(w)){
//                     //                 neg_neighbors2delete_weight.emplace(w, 0.0);
//                     //                 neg_neighbors2delete_count.emplace(w, 0);
//                     //             }
//                     //             neg_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                     //             neg_neighbors2delete_count[w] += 1;
//                     //         }
//                     //     }
//                     // }
//                     // for(auto &v_w : pos_neighbors2delete_weight){
//                     //     int v = v_w.first;
//                     //     double ori_weight = pos_edged_vertices[v] * graph[v].size();
//                     //     pos_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - pos_neighbors2delete_count[v]);
//                     // }
//                     // for(auto &v_w : neg_neighbors2delete_weight){
//                     //     int v = v_w.first;
//                     //     double ori_weight = neg_edged_vertices[v] * graph[v].size();
//                     //     neg_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - neg_neighbors2delete_count[v]);

//                     // }

//                     for (auto &v_s : pos_edged_vertices)
//                     {
//                         int v = v_s.first;
//                         double ori_w = pos_edged_vertices[v] * cc.size();
//                         for (int w : deleted_vertices)
//                         {
//                             Weight a = sim_graph[w][v] + 1;
//                             // ori_w -= log(a.numerator) - log(a.denominator);
//                             ori_w -= a.numerator / a.denominator;
//                         }
//                         // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
//                         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                     }
//                     for (auto &v_s : neg_edged_vertices)
//                     {
//                         int v = v_s.first;
//                         double ori_w = neg_edged_vertices[v] * cc.size();
//                         for (int w : deleted_vertices)
//                         {
//                             Weight a = sim_graph[w][v] + 1;
//                             ori_w -= a.numerator / a.denominator;
//                         }
//                         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                     }

//                     unordered_set<int> new_vertex;
//                     for (int u : deleted_vertices)
//                     {
//                         for (auto &v_sup : graph[u])
//                         {
//                             int v = v_sup.first;
//                             // graph[v].erase(u); // remove u from v's neighbors

//                             // check whether has new edged vertex
//                             if (!deleted_vertices.count(v) && !new_vertex.count(v) && !pos_edged_vertices.count(v) && !neg_edged_vertices.count(v))
//                             {
//                                 new_vertex.insert(v);
//                                 // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                             }
//                         }
//                     }
//                     // update graph
//                     for (int u : deleted_vertices)
//                     {
//                         unordered_map<int, int> neighbors = graph[u];
//                         vector<Edge> changed_edges;
//                         for (auto &v_sup : neighbors)
//                         {
//                             int v = v_sup.first;
//                             graph[v].erase(u);
//                             graph[u].erase(v);
//                             for (auto &it1 : graph[v])
//                             {
//                                 int w = it1.first;
//                                 if (graph[u].count(w))
//                                 {
//                                     --graph[v][w];
//                                     --graph[w][v];
//                                     changed_edges.push_back(pair<int, int>(v, w));
//                                 }
//                             }
//                         }

//                         cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
//                         graph.erase(u);
//                         VertexInfo removeV;
//                         removeV.vertex = u;
//                         removeV.neighbors = move(neighbors);
//                         removeV.changed_edges = move(changed_edges);
//                         Removed_vertices.push_back(removeV);
//                     }
//                     // printGraphAsEdgeList(graph, cc);
//                     // for (auto &u_v : deleted_edges)
//                     // {
//                     //     int u = u_v.first;
//                     //     for (auto &v_sup : u_v.second)
//                     //     {
//                     //         int v = v_sup.first;
//                     //         graph[u].erase(v);
//                     //     }
//                     // }

//                     for (auto &v : new_vertex)
//                     {
//                         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                     }

//                     // for (auto &u_v : changed_edges)
//                     // {
//                     //     int u = u_v.first;
//                     //     int v = u_v.second;
//                     //     graph[u][v]--;
//                     //     graph[v][u]--;
//                     // }

//                     break;
//                 }
//             }
//             if (u_it == vec.end())
//             {
//                 // vector<int> targets = {3, 25, 71, 168};
//                 // if (containsAll(cc, targets))
//                 // {
//                 // printGraphAsEdgeList(graph, cc);
//                 // }
//                 // Weight score = WEIGHT_MAX;
//                 // for (int i = 0; i < cc.size(); i++)
//                 // {
//                 //     int u = cc[i];

//                 //     for (int j = i + 1; j < cc.size(); j++)
//                 //     {
//                 //         int v = cc[j];
//                 //         if(sim_graph[u][v] == Weight(1, 204)){
//                 //             int a = 1;
//                 //         }
//                 //         score = min(score, sim_graph[u][v]);
//                 //     }
//                 // }
//                 // break;
//                 /*clear pos/neg set, find the non-cut vertex,
//                 select the vertex don't destroy the opposite sides,
//                 find the new pos/neg vertex from the deleted vertex*/
//                 // set<int> cv;
//                 // BCC(graph, cv);

//                 // connect all pos vertex and expand to k-truss
//                 // find the minimum degree neg vertex and remove it, check the remaining subgraph whether is connected
//                 vector<pair<int, int>> neg_vertices;
//                 for (auto &v : cc)
//                 {
//                     if (!seq2att[v].first)
//                     {
//                         neg_vertices.push_back(make_pair(v, graph[v].size()));
//                     }
//                 }
//                 sort(neg_vertices.begin(), neg_vertices.end(),
//                      [](const pair<int, int> &a, const pair<int, int> &b)
//                      {
//                          return a.second < b.second;
//                      });
//                 VertexInfo vertexInfo;
//                 auto u_it = neg_vertices.begin();
//                 // select a neg vertex
//                 for (; u_it != neg_vertices.end();)
//                 {
//                     int u = u_it->first;
//                     set<int> del_vertex;

//                     vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, del_vertex);
//                     del_vertex.insert(u);
//                     bool connected = false;
//                     for (int v : cc)
//                     {
//                         if (!del_vertex.count(v))
//                         {
//                             unordered_set<int> visited;
//                             DFS(graph, v, visited);

//                             if (visited.size() == cc.size() - del_vertex.size())
//                             {
//                                 connected = true;
//                             }
//                             break;
//                         }
//                     }

//                     if (!connected)
//                     {
//                         for (int v : del_vertex)
//                         {
//                             if (v != u)
//                             {
//                                 graph.emplace(v, unordered_map<int, int>());
//                             }
//                         }
//                         recoverTruss(vertexInfo, graph);
//                         u_it++;
//                         continue;
//                     }
//                     else
//                     {
//                         for (auto &v : del_vertex)
//                         {
//                             cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
//                         }
//                         int neg_num = 0;
//                         for (int v : del_vertex)
//                         {
//                             if (!seq2att[v].first)
//                             {
//                                 neg_num++;
//                             }
//                         }
//                         Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - del_vertex.size());
//                         if (new_ns > ns)
//                         {
//                             for (int v : del_vertex)
//                             {
//                                 if (v != u)
//                                 {
//                                     graph.emplace(v, unordered_map<int, int>());
//                                 }
//                             }
//                             recoverTruss(vertexInfo, graph);
//                             u_it++;
//                             continue;
//                         }
//                         else
//                         {
//                             ns = new_ns;
//                             // update pos/neg edged vertices
//                             pos_edged_vertices.clear();
//                             neg_edged_vertices.clear();

//                             unordered_set<int> new_vertex;
//                             for (auto &e : vertexInfo.sup_less_edges)
//                             {
//                                 int u = e.first.first;
//                                 int v = e.first.second;
//                                 if (!new_vertex.count(v) && !del_vertex.count(v))
//                                 {
//                                     new_vertex.insert(v);
//                                 }
//                                 if (!new_vertex.count(u) && !del_vertex.count(u))
//                                 {
//                                     new_vertex.insert(u);
//                                 }
//                             }
//                             // for (int u : del_vertex)
//                             // {
//                             //     for (auto &v_sup : graph[u])
//                             //     {
//                             //         int v = v_sup.first;
//                             //         // graph[v].erase(u); // remove u from v's neighbors

//                             //         // check whether has new edged vertex
//                             //         if (!del_vertex.count(v) && !new_vertex.count(v))
//                             //         {
//                             //             new_vertex.insert(v);
//                             //             // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                             //         }
//                             //     }
//                             // }
//                             for (auto &v : new_vertex)
//                             {
//                                 isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                             }

//                             break;
//                         }
//                     }
//                 }
//                 if (u_it == neg_vertices.end())
//                 {
//                     break;
//                 }
//             }
//         }
//         // // remove the neg vertices
//         // if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//         //                            seq2att, Removed_vertices,
//         //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//         // {
//         //     break;
//         //     // if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//         //     //                            seq2att, Removed_vertices,
//         //     //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//         //     // {
//         //     //     // has_vertex =  false;
//         //     //     break;
//         //     // }
//         // }
//     }

//     if (ns >= thd1 && ns <= thd2)
//     {
//         Weight score = WEIGHT_MAX;
//         for (int i = 0; i < cc.size(); i++)
//         {
//             int u = cc[i];

//             for (int j = i + 1; j < cc.size(); j++)
//             {
//                 int v = cc[j];

//                 score = min(score, sim_graph[u][v]);
//             }
//         }
//         // cout << "score: " << score.numerator << " / " << score.denominator << endl;
//         // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
//         if (score > bottom_score)
//         {
//             bottom_score = score;
//             if (bottom_score == Weight(1, 526))
//             // if (bottom_score == Weight(1, 192))
//             // if (bottom_score == Weight(1, 6336))
//             {
//                 // printGraphAsEdgeList(graph, cc);
//                 int a = 1;
//             }
//             cout << "bottom_score: " << bottom_score.numerator << "/" << bottom_score.denominator << endl;
//             printGraphAsEdgeList(graph, cc);
//             // remove vertex with score less than bottom_score and maintain the truss
//             vector<VertexInfo> changed_edges;
//             vector<set<int>> del_vertices;
//             greedy_remove(changed_edges, del_vertices, cc, graph, sim_graph, bottom_score, k);
//             // find connected components
//             vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
//             // calculate all added vertex-pair score
//             // set<int> added_vertices;
//             for (auto &cc : connectedComponents)
//             {
//                 if (cc.size() <= 1)
//                     continue;
//                 // repeat the process
//                 // printGraphAsEdgeList(graph, cc);
//                 cout << "delete edges" << endl;
//                 shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//             }
//             for (auto &it : del_vertices)
//             {
//                 for (int v : it)
//                 {
//                     graph.emplace(v, unordered_map<int, int>());
//                 }
//             }

//             while (!changed_edges.empty())
//             {
//                 auto &it = changed_edges.back();

//                 recoverTruss(it, graph);
//                 changed_edges.pop_back();
//             }
//         }
//     }

//     // recovery
//     // if (score_id == score_edges.size())
//     //     break;
//     while (!Removed_vertices.empty())
//     {
//         VertexInfo v = Removed_vertices.back();
//         RecoverEdges(v.vertex, v.neighbors, v.changed_edges, graph);
//         Removed_vertices.pop_back();
//     }
// }
// 5-27
void shrink(std::unordered_map<int, std::unordered_map<int, int>> &graph, vector<int> &cc, Weight &bottom_score,
            map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
            unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    Weight ns = NSize(cc, seq2att);
    if (cc.size() == k)
    {
        if (ns >= thd1 && ns <= thd2)
        {
            Weight score = WEIGHT_MAX;
            for (int i = 0; i < cc.size(); i++)
            {
                int u = cc[i];

                for (int j = i + 1; j < cc.size(); j++)
                {
                    int v = cc[j];

                    score = min(score, sim_graph[u][v]);
                }
            }
            bottom_score = score;
            // printGraphAsEdgeList(graph, cc);
        }
        return;
    }
    if (thd2 < 1 && ns == Weight(1, 1))
    {
        return;
    }
    if (thd1 > 0 && ns == Weight(0, 1))
    {
        return;
    }
    vector<VertexInfo> Removed_vertices;
    // printGraphAsEdgeList(graph, cc);
    // unordered_set<int> pos_edged_vertices;
    // unordered_set<int> neg_edged_vertices;
    unordered_map<int, double> pos_edged_vertices;
    unordered_map<int, double> neg_edged_vertices;
    for (int v : cc)
    {
        // if(!graph.count(src)) continue;
        // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
        double weight = calcul_weight(graph, cc, v, sim_graph);
        // Weight weight = calcul_weight(graph, cc, v, sim_graph);
        if (seq2att[v].first)
        {
            pos_edged_vertices.emplace(v, weight);
        }
        else
        {
            // neg_edged_vertices.emplace(v, sup);
            neg_edged_vertices.emplace(v, weight);
        }
    }
    // double ns = NSize(edge_truss, seq2att);
    // double ns = 0.0;
    // for (int u : cc)
    // {
    //     if (!seq2att[u].first)
    //     {
    //         ns++;
    //     }
    // }
    // ns /= cc.size();
    while (cc.size() > k && (ns < thd1 || ns > thd2))
    {
        // bool has_vertex = true;

        if (ns < thd1)
        {
            unordered_map<int, double> new_pos_edged_vertices;
            for (auto &u_w : pos_edged_vertices)
            {
                int u = u_w.first;
                int neg_neighbor_size = 0;
                for (auto &u_v : graph[u])
                {
                    int v = u_v.first;
                    if (!seq2att[v].first)
                    {
                        neg_neighbor_size++;
                    }
                }
                new_pos_edged_vertices.emplace(u, u_w.second + neg_neighbor_size);
            }
            // vector<pair<int, double>> vec(neg_edged_vertices.begin(), neg_edged_vertices.end());
            vector<pair<int, double>> vec(new_pos_edged_vertices.begin(), new_pos_edged_vertices.end());
            sort(vec.begin(), vec.end(), compareByValue);
            VertexInfo vertexInfo;
            auto u_it = vec.begin();
            for (; u_it != vec.end();)
            {
                int u = u_it->first;
                set<int> del_vertex;

                vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, del_vertex);

                del_vertex.insert(u);
                bool connected = false;
                for (int v : cc)
                {
                    if (!del_vertex.count(v))
                    {
                        unordered_set<int> visited;
                        DFS(graph, v, visited);

                        if (visited.size() == cc.size() - del_vertex.size())
                        {
                            connected = true;
                        }
                        break;
                    }
                }

                if (!connected)
                {
                    for (int v : del_vertex)
                    {
                        if (v != u)
                        {
                            graph.emplace(v, unordered_map<int, int>());
                        }
                    }
                    recoverTruss(vertexInfo, graph);
                    // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                    u_it++;
                    continue;
                }
                else
                {
                    int neg_num = 0;
                    for (int v : del_vertex)
                    {
                        if (!seq2att[v].first)
                        {
                            neg_num++;
                        }
                    }
                    Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - del_vertex.size());
                    if (new_ns < ns)
                    {
                        for (int v : del_vertex)
                        {
                            if (v != u)
                            {
                                graph.emplace(v, unordered_map<int, int>());
                            }
                        }
                        recoverTruss(vertexInfo, graph);
                        // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                        u_it++;
                        continue;
                    }
                    else
                    {
                        ns = new_ns;
                        Removed_vertices.push_back(vertexInfo);
                        for (auto &v : del_vertex)
                        {
                            cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
                            if (pos_edged_vertices.count(v))
                            {
                                pos_edged_vertices.erase(v);
                            }
                            else if (neg_edged_vertices.count(v))
                            {
                                neg_edged_vertices.erase(v);
                            }
                        }
                        // printGraphAsEdgeList(graph, cc);
                        // update pos/neg edged vertices

                        unordered_set<int> new_vertex;
                        for (auto &e : vertexInfo.sup_less_edges)
                        {
                            int u = e.first.first;
                            int v = e.first.second;
                            if (!new_vertex.count(v) && !del_vertex.count(v))
                            {
                                new_vertex.insert(v);
                            }
                            if (!new_vertex.count(u) && !del_vertex.count(u))
                            {
                                new_vertex.insert(u);
                            }
                        }

                        for (auto &v : new_vertex)
                        {
                            if (pos_edged_vertices.count(v))
                            {
                                double ori_w = pos_edged_vertices[v] * cc.size();
                                for (int w : del_vertex)
                                {
                                    Weight a = sim_graph[w][v] + 1;
                                    // ori_w -= log(a.numerator) - log(a.denominator);
                                    ori_w -= a.numerator / a.denominator;
                                }
                                // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
                                pos_edged_vertices[v] = ori_w / (cc.size() - del_vertex.size());
                            }
                            else if (neg_edged_vertices.count(v))
                            {
                                double ori_w = neg_edged_vertices[v] * cc.size();
                                for (int w : del_vertex)
                                {
                                    Weight a = sim_graph[w][v] + 1;
                                    // ori_w -= log(a.numerator) - log(a.denominator);
                                    ori_w -= a.numerator / a.denominator;
                                }
                                // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
                                neg_edged_vertices[v] = ori_w / (cc.size() - del_vertex.size());
                            }
                        }

                        break;
                    }
                }
            }
            if (u_it == vec.end())
            {
                break;
            }
        }
        else if (ns > thd2)
        {
            unordered_map<int, double> new_neg_edged_vertices;
            for (auto &u_w : neg_edged_vertices)
            {
                int u = u_w.first;
                int pos_neighbor_size = 0;
                for (auto &u_v : graph[u])
                {
                    int v = u_v.first;
                    if (seq2att[v].first)
                    {
                        pos_neighbor_size++;
                    }
                }
                new_neg_edged_vertices.emplace(u, u_w.second + pos_neighbor_size);
            }
            // vector<pair<int, double>> vec(neg_edged_vertices.begin(), neg_edged_vertices.end());
            vector<pair<int, double>> vec(new_neg_edged_vertices.begin(), new_neg_edged_vertices.end());
            sort(vec.begin(), vec.end(), compareByValue);
            VertexInfo vertexInfo;
            auto u_it = vec.begin();
            for (; u_it != vec.end();)
            {
                int u = u_it->first;
                set<int> del_vertex;

                vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, del_vertex);

                del_vertex.insert(u);
                bool connected = false;
                for (int v : cc)
                {
                    if (!del_vertex.count(v))
                    {
                        unordered_set<int> visited;
                        DFS(graph, v, visited);

                        if (visited.size() == cc.size() - del_vertex.size())
                        {
                            connected = true;
                        }
                        break;
                    }
                }

                if (!connected)
                {
                    for (int v : del_vertex)
                    {
                        if (v != u)
                        {
                            graph.emplace(v, unordered_map<int, int>());
                        }
                    }
                    recoverTruss(vertexInfo, graph);
                    // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                    u_it++;
                    continue;
                }
                else
                {
                    int neg_num = 0;
                    for (int v : del_vertex)
                    {
                        if (!seq2att[v].first)
                        {
                            neg_num++;
                        }
                    }
                    Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - del_vertex.size());
                    if (new_ns > ns)
                    {
                        for (int v : del_vertex)
                        {
                            if (v != u)
                            {
                                graph.emplace(v, unordered_map<int, int>());
                            }
                        }
                        recoverTruss(vertexInfo, graph);
                        // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                        u_it++;
                        continue;
                    }
                    else
                    {
                        ns = new_ns;
                        Removed_vertices.push_back(vertexInfo);
                        for (auto &v : del_vertex)
                        {
                            cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
                            if (pos_edged_vertices.count(v))
                            {
                                pos_edged_vertices.erase(v);
                            }
                            else if (neg_edged_vertices.count(v))
                            {
                                neg_edged_vertices.erase(v);
                            }
                        }
                        // printGraphAsEdgeList(graph, cc);
                        // update pos/neg edged vertices

                        unordered_set<int> new_vertex;
                        for (auto &e : vertexInfo.sup_less_edges)
                        {
                            int u = e.first.first;
                            int v = e.first.second;
                            if (!new_vertex.count(v) && !del_vertex.count(v))
                            {
                                new_vertex.insert(v);
                            }
                            if (!new_vertex.count(u) && !del_vertex.count(u))
                            {
                                new_vertex.insert(u);
                            }
                        }

                        for (auto &v : new_vertex)
                        {
                            if (pos_edged_vertices.count(v))
                            {
                                double ori_w = pos_edged_vertices[v] * cc.size();
                                for (int w : del_vertex)
                                {
                                    Weight a = sim_graph[w][v] + 1;
                                    // ori_w -= log(a.numerator) - log(a.denominator);
                                    ori_w -= a.numerator / a.denominator;
                                }
                                // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
                                pos_edged_vertices[v] = ori_w / (cc.size() - del_vertex.size());
                            }
                            else if (neg_edged_vertices.count(v))
                            {
                                double ori_w = neg_edged_vertices[v] * cc.size();
                                for (int w : del_vertex)
                                {
                                    Weight a = sim_graph[w][v] + 1;
                                    // ori_w -= log(a.numerator) - log(a.denominator);
                                    ori_w -= a.numerator / a.denominator;
                                }
                                // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
                                neg_edged_vertices[v] = ori_w / (cc.size() - del_vertex.size());
                            }
                        }

                        break;
                    }
                }
            }
            if (u_it == vec.end())
            {
                break;
            }
        }
    }

    if (ns >= thd1 && ns <= thd2)
    {
        Weight score = WEIGHT_MAX;
        for (int i = 0; i < cc.size(); i++)
        {
            int u = cc[i];

            for (int j = i + 1; j < cc.size(); j++)
            {
                int v = cc[j];

                score = min(score, sim_graph[u][v]);
            }
        }
        // cout << "score: " << score.numerator << " / " << score.denominator << endl;
        // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
        if (score > bottom_score)
        {
            bottom_score = score;
            if (bottom_score == Weight(1, 526))
            // if (bottom_score == Weight(1, 192))
            // if (bottom_score == Weight(1, 6336))
            {
                // printGraphAsEdgeList(graph, cc);
                int a = 1;
            }
            cout << "bottom_score: " << bottom_score.numerator << "/" << bottom_score.denominator << endl;
            // printGraphAsEdgeList(graph, cc);
            // remove vertex with score less than bottom_score and maintain the truss
            vector<VertexInfo> changed_edges;
            vector<set<int>> del_vertices;
            greedy_remove(changed_edges, del_vertices, cc, graph, sim_graph, bottom_score, k);
            // find connected components
            vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
            // calculate all added vertex-pair score
            // set<int> added_vertices;
            for (auto &cc : connectedComponents)
            {
                if (cc.size() <= 1)
                    continue;
                // repeat the process
                // printGraphAsEdgeList(graph, cc);
                cout << "delete edges" << endl;
                shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
            }
            for (auto &it : del_vertices)
            {
                for (int v : it)
                {
                    graph.emplace(v, unordered_map<int, int>());
                }
            }

            while (!changed_edges.empty())
            {
                auto &it = changed_edges.back();

                recoverTruss(it, graph);
                changed_edges.pop_back();
            }
        }
    }

    // recovery
    // if (score_id == score_edges.size())
    //     break;
    while (!Removed_vertices.empty())
    {
        VertexInfo v = Removed_vertices.back();
        // RecoverEdges(v.vertex, v.neighbors, v.changed_edges, graph);
        recoverTruss(v, graph);
        Removed_vertices.pop_back();
    }
}
// bool special_vertex(int u, unordered_map<int, unordered_map<int, int>> &graph, int k,
//                     set<int> &deleted_vertices,
//                     unordered_map<int, unordered_map<int, int>> &deleted_edges)
// {
//     // get 1-hop subgraph
//     unordered_map<int, unordered_map<int, int>> subgraph;
//     const auto &neighbors = graph[u];
//     for (const auto &[v, weight_v] : neighbors)
//     {
//         subgraph[v] = {};
//         // subgraph[u][v] = weight_v; // 添加 u -> v 边
//         for (const auto &[w, weight_w] : graph.at(v))
//         {
//             // 只保留邻居之间的连接
//             if (neighbors.count(w))
//             {
//                 subgraph[v][w] = weight_w;
//             }
//         }
//     }

//     unordered_set<int> visited;
//     // 选择 neighbors 中的一个节点作为起始节点
//     int startNode = neighbors.begin()->first;

//     // 从起始节点开始进行 DFS 遍历
//     DFS(subgraph, startNode, visited);

//     // 检查 neighbors 中的所有节点是否都被访问过
//     if (visited.size() != neighbors.size())
//     {
//         return false;
//     }
//     unordered_set<int> deleted_vertices_cand;
//     for (auto &u_v_w : subgraph)
//     {
//         int u = u_v_w.first;
//         if (subgraph[u].size() == graph[u].size() - 1)
//         {
//             deleted_vertices_cand.insert(u);
//         }
//     }
//     // unordered_map<int, unordered_map<int, int>> deleted_edges;
//     for (auto it = neighbors.begin(); it != neighbors.end(); it++)
//     {
//         int v = it->first;

//         for (auto &it1 : subgraph[v])
//         {

//             int w = it1.first;
//             if (v > w)
//                 continue;
//             // int new_sup = -1;
//             // update support
//             // 需要证明如果边缘点候选存在临边是被删除检查点影响的边，那么检查点的删除不影响这个点是否还是边缘点。具体的，他们的关系是共同构成一个团
//             int sup = subgraph[v][w];
//             if (sup == k - 2)
//             {
//                 if (!deleted_vertices_cand.count(w) && !deleted_vertices_cand.count(v))
//                 {
//                     return false;
//                 }
//                 if (!deleted_edges.count(v))
//                 {
//                     deleted_edges.emplace(v, unordered_map<int, int>());
//                     // deleted_edges.emplace(v, unordered_set<int>());
//                 }
//                 deleted_edges[v].emplace(w, sup);
//                 // deleted_edges[v].insert(w);
//                 if (!deleted_edges.count(w))
//                 {
//                     deleted_edges.emplace(w, unordered_map<int, int>());
//                     // deleted_edges.emplace(w, unordered_set<int>());
//                 }
//                 deleted_edges[w].emplace(v, sup);
//                 // deleted_edges[w].insert(v);
//             }
//             // else
//             // {
//             //     changed_edges.push_back(Edge(v, w));
//             // }
//         }
//     }
//     // unordered_set<int> deleted_vertices;
//     for (auto &u_v : deleted_edges)
//     {
//         int u = u_v.first;
//         // if (subgraph[u].size() != graph[u].size() - 1)
//         // {
//         //     deleted_vertices.clear();
//         //     // deleted_edges.clear();
//         //     return;
//         // }
//         if (deleted_vertices_cand.count(u) && deleted_edges[u].size() == subgraph[u].size())
//         {
//             deleted_vertices.insert(u);
//         }
//     }

//     for (auto &u : deleted_vertices)
//     {

//         for (auto &v_s : subgraph[u])
//         {
//             int v = v_s.first;
//             subgraph[v].erase(u);
//         }
//         subgraph.erase(u);
//     }
//     if (subgraph.empty())
//     {
//         deleted_vertices.clear();
//         deleted_edges.clear();
//         return false;
//     }
//     // unordered_set<int> visited;
//     visited.clear();
//     int w = subgraph.begin()->first;
//     DFS(subgraph, w, visited);
//     if (visited.size() == subgraph.size())
//     {
//         deleted_vertices.insert(u);
//         return true;
//         // deleted_edges.emplace(u, unordered_map<int, int>());

//         // for (auto &v_sup : graph[u])
//         // {
//         //     int v = v_sup.first;
//         //     int sup = v_sup.second;
//         //     deleted_edges[u].emplace(v, sup);
//         //     if (!deleted_edges.count(v))
//         //     {
//         //         deleted_edges.emplace(v, unordered_map<int, int>());
//         //         // deleted_edges.emplace(v, unordered_set<int>());
//         //     }
//         //     deleted_edges[v].emplace(w, sup);
//         //     // deleted_edges[v].insert(w);
//         // }
//     }
//     else
//     {
//         deleted_vertices.clear();
//         deleted_edges.clear();
//         return false;
//         // changed_edges.clear();
//     }
//     // if (visited.size() == subgraph.size())
//     // {
//     //     // deleted_vertices.insert(u);
//     //     // return move(deleted_vertices);
//     // }
//     // else
//     // {
//     //     return {};
//     // }
// }
// bool special_vertex(int u, unordered_map<int, unordered_map<int, int>> &graph, int k,
//                     set<int> &deleted_vertices)
// {
//     // get 1-hop subgraph
//     unordered_map<int, unordered_map<int, int>> subgraph;
//     const auto &neighbors = graph[u];
//     for (const auto &[v, weight_v] : neighbors)
//     {
//         subgraph[v] = {};
//         // subgraph[u][v] = weight_v; // 添加 u -> v 边
//         for (const auto &[w, weight_w] : graph.at(v))
//         {
//             // 只保留邻居之间的连接
//             if (neighbors.count(w))
//             {
//                 subgraph[v][w] = weight_w;
//             }
//         }
//     }

//     unordered_set<int> visited;
//     // 选择 neighbors 中的一个节点作为起始节点
//     int startNode = neighbors.begin()->first;

//     // 从起始节点开始进行 DFS 遍历
//     DFS(subgraph, startNode, visited);

//     // 检查 neighbors 中的所有节点是否都被访问过
//     if (visited.size() != neighbors.size())
//     {
//         return false;
//     }
//     unordered_set<int> deleted_vertices_cand;
//     for (auto &u_v_w : subgraph)
//     {
//         int u = u_v_w.first;
//         if (subgraph[u].size() == graph[u].size() - 1)
//         {
//             deleted_vertices_cand.insert(u);
//         }
//     }
//     // unordered_map<int, unordered_map<int, int>> deleted_edges;
//     vector<Edge> deleted_edges;
//     for (auto it = neighbors.begin(); it != neighbors.end(); it++)
//     {
//         int v = it->first;

//         for (auto &it1 : subgraph[v])
//         {

//             int w = it1.first;
//             if (v > w)
//                 continue;
//             // int new_sup = -1;
//             // update support
//             // 需要证明如果边缘点候选存在临边是被删除检查点影响的边，那么检查点的删除不影响这个点是否还是边缘点。具体的，他们的关系是共同构成一个团
//             int &sup = subgraph[v][w];
//             if (sup == k - 2)
//             {
//                 if (!deleted_vertices_cand.count(w) && !deleted_vertices_cand.count(v))
//                 {
//                     return false;
//                 }
//                 deleted_edges.push_back(Edge(v, w));
//             }
//             sup--;
//             subgraph[w][v] = sup;
//             // else
//             // {
//             //     changed_edges.push_back(Edge(v, w));
//             // }
//         }
//     }
//     // vector<Edge> deleted_edges;
//     size_t i = 0;
//     for (; i < deleted_edges.size(); i++)
//     {

//         Edge e = deleted_edges[i];
//         int u = e.first, v = e.second;

//         // 确保u的度数小于等于v的度数
//         if (subgraph[u].size() > subgraph[v].size())
//             swap(u, v);

//         // 遍历u的邻居
//         for (const auto &nbr_entry : subgraph[u])
//         {
//             int nbr = nbr_entry.first;
//             if (!subgraph[v].count(nbr))
//             {
//                 continue;
//             }

//             Edge nbr_v = GetEdge(nbr, v);
//             Edge nbr_u = GetEdge(nbr, u);

//             // 更新nbr_v边的支持度
//             if (subgraph[nbr][v] != -1)
//             {
//                 int &support = subgraph[nbr][v];
//                 if (support == k - 2)
//                 {
//                     deleted_edges.push_back(nbr_v);
//                 }
//                 support--;
//                 subgraph[v][nbr] = support;
//             }

//             // 更新nbr_u边的支持度
//             if (subgraph[nbr][u] != -1)
//             {
//                 int &support = subgraph[nbr][u];
//                 if (support == k - 2)
//                 {
//                     deleted_edges.push_back(nbr_u);
//                 }
//                 support--;
//                 subgraph[u][nbr] = support;
//             }
//         }

//         // 从图中删除这条边
//         subgraph[u].erase(v);
//         subgraph[v].erase(u);
//     }
//     // get isloated vertices
//     for (auto &u_v : subgraph)
//     {
//         int v = u_v.first;
//         if (u_v.second.size() == 0)
//         {
//             if (deleted_vertices_cand.count(u))
//             {
//                 deleted_vertices.insert(u_v.first);
//                 subgraph.erase(v);
//             }
//             else
//             {
//                 deleted_vertices.clear();
//                 return false;
//             }
//         }
//     }

//     if (subgraph.empty())
//     {
//         deleted_vertices.clear();
//         deleted_edges.clear();
//         return false;
//     }
//     // unordered_set<int> visited;
//     visited.clear();
//     int w = subgraph.begin()->first;
//     DFS(subgraph, w, visited);
//     if (visited.size() == subgraph.size())
//     {
//         deleted_vertices.insert(u);
//         return true;
//         // deleted_edges.emplace(u, unordered_map<int, int>());

//         // for (auto &v_sup : graph[u])
//         // {
//         //     int v = v_sup.first;
//         //     int sup = v_sup.second;
//         //     deleted_edges[u].emplace(v, sup);
//         //     if (!deleted_edges.count(v))
//         //     {
//         //         deleted_edges.emplace(v, unordered_map<int, int>());
//         //         // deleted_edges.emplace(v, unordered_set<int>());
//         //     }
//         //     deleted_edges[v].emplace(w, sup);
//         //     // deleted_edges[v].insert(w);
//         // }
//     }
//     else
//     {
//         deleted_vertices.clear();
//         return false;
//         // changed_edges.clear();
//     }
//     // if (visited.size() == subgraph.size())
//     // {
//     //     // deleted_vertices.insert(u);
//     //     // return move(deleted_vertices);
//     // }
//     // else
//     // {
//     //     return {};
//     // }
// }
void update(unordered_map<int, double> &pos_vertices, unordered_map<int, double> &neg_vertices,
            set<int> &deleted_vertices, vector<int> &cc, unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    for (int u : deleted_vertices)
    {
        if (pos_vertices.count(u))
        {
            pos_vertices.erase(u);
        }
        else if (neg_vertices.count(u))
        {
            neg_vertices.erase(u);
        }
    }

    for (auto &v_s : pos_vertices)
    {
        int v = v_s.first;
        double ori_w = pos_vertices[v] * (cc.size() - 1);
        for (int w : deleted_vertices)
        {
            // Weight a = sim_graph[w][v] + 1;
            Weight a = sim_graph[w][v];
            // ori_w -= log(a.numerator) - log(a.denominator);
            ori_w -= double(a.numerator) / a.denominator;
        }
        // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
        v_s.second = ori_w / (cc.size() - deleted_vertices.size() - 1);
    }
    for (auto &v_s : neg_vertices)
    {
        int v = v_s.first;
        double ori_w = neg_vertices[v] * cc.size();
        for (int w : deleted_vertices)
        {
            Weight a = sim_graph[w][v];
            ori_w -= double(a.numerator) / a.denominator;
        }
        v_s.second = ori_w / (cc.size() - deleted_vertices.size() - 1);
    }
}
void vertex_weight_aug(bool sign, unordered_map<int, unordered_map<int, int>> &graph,
                       unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
                       unordered_map<int, double> &vertices, unordered_map<int, double> &new_vertices)
{
    for (auto &u_w : vertices)
    {
        int u = u_w.first;
        int neg_neighbor_size = 0;
        for (auto &u_v : graph[u])
        {
            int v = u_v.first;
            if (sign & seq2att[v].first)
            {
                neg_neighbor_size++;
            }
        }
        new_vertices.emplace(u, u_w.second + neg_neighbor_size);
    }
}
Weight ns_after_delete(Weight &ns, set<int> &del_vertex, vector<int> &cc,
                       unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    int neg_num = 0;
    for (int v : del_vertex)
    {
        if (!seq2att[v].first)
        {
            neg_num++;
        }
    }
    Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - del_vertex.size());
    return move(new_ns);
}
bool edge_vertex(int u, unordered_map<int, unordered_map<int, int>> &graph, int k, VertexInfo &vertexInfo)
{
    vertexInfo.vertex = u;
    // delete vertex
    vertexInfo.neighbors = graph.at(u);
    // get 1-hop subgraph
    unordered_map<int, unordered_map<int, int>> subgraph;
    // const auto &neighbors = graph[u];
    for (const auto &[v, weight_v] : vertexInfo.neighbors)
    {
        subgraph[v] = {};
        // subgraph[u][v] = weight_v; // 添加 u -> v 边
        for (const auto &[w, weight_w] : graph.at(v))
        {
            // 只保留邻居之间的连接
            if (vertexInfo.neighbors.count(w))
            {
                subgraph[v][w] = weight_w;
            }
        }
    }

    // delete u from kct

    unordered_set<int> visited;
    // 选择 neighbors 中的一个节点作为起始节点
    int startNode = vertexInfo.neighbors.begin()->first;

    // 从起始节点开始进行 DFS 遍历
    DFS(subgraph, startNode, visited);

    // 检查 neighbors 中的所有节点是否都被访问过
    if (visited.size() != vertexInfo.neighbors.size())
    {
        return false;
    }
    for (auto &v_w_s : subgraph)
    {
        int v = v_w_s.first;
        for (auto &w_s : subgraph[v])
        {
            int w = w_s.first;
            if (w_s.second == k - 2)
            {
                // for(auto &e : vertexInfo.changed_edges){
                //     int v = e.first;
                //     int w = e.second;
                //     graph[v][w]++;
                // }
                // vertexInfo.neighbors.clear();
                // vertexInfo.changed_edges.clear();
                return false;
            }
            // w_s.second--;
            // graph[v][w]--;
            // vertexInfo.changed_edges.push_back(Edge(v, w));
        }
    }
    // printGraphAsEdgeList(graph);
    for (auto &v_w_s : subgraph)
    {
        int v = v_w_s.first;
        for (auto &w_s : subgraph[v])
        {
            int w = w_s.first;
            graph[v][w]--;
            vertexInfo.changed_edges.push_back(Edge(v, w));
        }
    }
    graph.erase(u);
    for (const auto &it : vertexInfo.neighbors)
    {
        if (graph.count(it.first))
        {
            graph[it.first].erase(u);
        }
    }
    return true;
}
void special_shrink(bool special, unordered_map<int, unordered_map<int, int>> &graph, vector<int> &cc, Weight &bottom_score,
                    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
                    unordered_map<int, unordered_map<int, Weight>> &sim_graph,
                    chrono::high_resolution_clock::time_point startTime)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - startTime;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return;
    }
    Weight ns = NSize(cc, seq2att);
    if (cc.size() == k)
    {
        if (ns >= thd1 && ns <= thd2)
        {
            Weight score = WEIGHT_MAX;
            for (int i = 0; i < cc.size(); i++)
            {
                int u = cc[i];

                for (int j = i + 1; j < cc.size(); j++)
                {
                    int v = cc[j];

                    score = min(score, sim_graph[u][v]);
                }
            }
            bottom_score = score;
            // printGraphAsEdgeList(graph, cc);
        }
        return;
    }
    if (thd2 < 1 && ns == Weight(1, 1))
    {
        return;
    }
    if (thd1 > 0 && ns == Weight(0, 1))
    {
        return;
    }
    vector<VertexInfo> Removed_vertices;
    // printGraphAsEdgeList(graph, cc);
    // unordered_set<int> pos_edged_vertices;
    // unordered_set<int> neg_edged_vertices;
    unordered_map<int, double> pos_vertices;
    unordered_map<int, double> neg_vertices;
    for (int v : cc)
    {
        // if(!graph.count(src)) continue;
        // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
        double weight = calcul_weight(graph, cc, v, sim_graph);
        // Weight weight = calcul_weight(graph, cc, v, sim_graph);
        if (seq2att[v].first)
        {
            pos_vertices.emplace(v, weight);
        }
        else
        {
            // neg_edged_vertices.emplace(v, sup);
            neg_vertices.emplace(v, weight);
        }
    }
    // double ns = NSize(edge_truss, seq2att);
    // double ns = 0.0;
    // for (int u : cc)
    // {
    //     if (!seq2att[u].first)
    //     {
    //         ns++;
    //     }
    // }
    // ns /= cc.size();
    while (cc.size() > k && (ns < thd1 || ns > thd2))
    {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - startTime;
        // cout << "Time used: " << duration.count() << "s" << endl;
        if (duration.count() > max_time)
        {
            return;
        }
        // bool has_vertex = true;

        if (ns < thd1)
        {
            unordered_map<int, double> new_pos_vertices;
            vertex_weight_aug(false, graph, seq2att, pos_vertices, new_pos_vertices);
            // vector<pair<int, double>> vec(neg_edged_vertices.begin(), neg_edged_vertices.end());
            vector<pair<int, double>> vec(new_pos_vertices.begin(), new_pos_vertices.end());
            sort(vec.begin(), vec.end(), compareByValue);

            auto u_it = vec.begin();
            for (; u_it != vec.end();)
            {
                VertexInfo vertexInfo;
                int u = u_it->first;
                set<int> deleted_vertices;
                // unordered_map<int, unordered_map<int, int>> deleted_edges;
                unordered_map<int, unordered_map<int, int>> deleted_edges;
                // if (special_vertex(u, graph, k, deleted_vertices, deleted_edges))
                if (Removed_vertices.size() == 9 && u == 22)
                {
                    int a = 1;
                }
                if (special && edge_vertex(u, graph, k, vertexInfo))
                {
                    ns = (ns * cc.size()) / (cc.size() - 1);

                    // update pos/neg edged vertices
                    deleted_vertices.insert(u);
                    // ns = ns_after_delete(ns, deleted_vertices, cc, seq2att);
                    update(pos_vertices, neg_vertices, deleted_vertices, cc, sim_graph);
                    // update graph
                    // set<int> del_vertex;
                    // deleted_vertices.clear();
                    // vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, deleted_vertices);
                    // deleted_vertices.insert(u);
                    Removed_vertices.push_back(vertexInfo);
                    cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
                    break;
                }
                else
                {
                    vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, deleted_vertices);
                    deleted_vertices.insert(u);
                    bool connected = false;
                    for (int v : cc)
                    {
                        if (!deleted_vertices.count(v))
                        {
                            unordered_set<int> visited;
                            DFS(graph, v, visited);

                            if (visited.size() == cc.size() - deleted_vertices.size())
                            {
                                connected = true;
                            }
                            break;
                        }
                    }

                    if (!connected)
                    {
                        // for (int v : deleted_vertices)
                        // {
                        //     if (v != u)
                        //     {
                        //         graph.emplace(v, unordered_map<int, int>());
                        //     }
                        // }
                        recoverTruss(vertexInfo, graph);
                        // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                        u_it++;
                        continue;
                    }
                    else
                    {
                        Weight new_ns = ns_after_delete(ns, deleted_vertices, cc, seq2att);
                        if (new_ns < ns)
                        {
                            // for (int v : deleted_vertices)
                            // {
                            //     if (v != u)
                            //     {
                            //         graph.emplace(v, unordered_map<int, int>());
                            //     }
                            // }
                            recoverTruss(vertexInfo, graph);
                            // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                            u_it++;
                            continue;
                        }
                        else
                        {
                            ns = new_ns;
                            update(pos_vertices, neg_vertices, deleted_vertices, cc, sim_graph);
                            Removed_vertices.push_back(vertexInfo);
                            for (auto &v : deleted_vertices)
                            {
                                cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
                            }

                            break;
                        }
                    }
                }
            }
            if (u_it == vec.end())
            {
                break;
            }
        }
        else if (ns > thd2)
        {
            unordered_map<int, double> new_neg_vertices;
            vertex_weight_aug(true, graph, seq2att, neg_vertices, new_neg_vertices);
            // for (auto &u_w : neg_vertices)
            // {
            //     int u = u_w.first;
            //     int pos_neighbor_size = 0;
            //     for (auto &u_v : graph[u])
            //     {
            //         int v = u_v.first;
            //         if (seq2att[v].first)
            //         {
            //             pos_neighbor_size++;
            //         }
            //     }
            //     new_neg_vertices.emplace(u, u_w.second + pos_neighbor_size);
            // }
            vector<pair<int, double>> vec(new_neg_vertices.begin(), new_neg_vertices.end());
            sort(vec.begin(), vec.end(), compareByValue);

            auto u_it = vec.begin();
            for (; u_it != vec.end();)
            {
                VertexInfo vertexInfo;
                int u = u_it->first;
                set<int> deleted_vertices;
                // unordered_map<int, unordered_map<int, int>> deleted_edges;
                unordered_map<int, unordered_map<int, int>> deleted_edges;
                // if (special_vertex(u, graph, k, deleted_vertices, deleted_edges))
                if (special && edge_vertex(u, graph, k, vertexInfo))
                {
                    // printGraphAsEdgeList(graph, cc);
                    ns = (ns * cc.size() - 1) / (cc.size() - 1);
                    // update pos/neg edged vertices
                    deleted_vertices.insert(u);
                    // ns = ns_after_delete(ns, deleted_vertices, cc, seq2att);
                    update(pos_vertices, neg_vertices, deleted_vertices, cc, sim_graph);
                    // update graph
                    // set<int> del_vertex;
                    // deleted_vertices.clear();
                    // vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, deleted_vertices);
                    // deleted_vertices.insert(u);
                    Removed_vertices.push_back(vertexInfo);
                    cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
                    break;
                }
                else
                {

                    vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, deleted_vertices);
                    deleted_vertices.insert(u);
                    bool connected = false;
                    for (int v : cc)
                    {
                        if (!deleted_vertices.count(v))
                        {
                            unordered_set<int> visited;
                            DFS(graph, v, visited);

                            if (visited.size() == cc.size() - deleted_vertices.size())
                            {
                                connected = true;
                            }
                            break;
                        }
                    }

                    if (!connected)
                    {
                        // for (int v : deleted_vertices)
                        // {
                        //     if (v != u)
                        //     {
                        //         graph.emplace(v, unordered_map<int, int>());
                        //     }
                        // }
                        recoverTruss(vertexInfo, graph);
                        // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                        u_it++;
                        continue;
                    }
                    else
                    {
                        Weight new_ns = ns_after_delete(ns, deleted_vertices, cc, seq2att);
                        if (new_ns > ns)
                        {
                            // for (int v : deleted_vertices)
                            // {
                            //     if (v != u)
                            //     {
                            //         graph.emplace(v, unordered_map<int, int>());
                            //     }
                            // }
                            recoverTruss(vertexInfo, graph);
                            // RecoverEdges(vertexInfo.vertex, vertexInfo.neighbors, vertexInfo.changed_edges, graph);
                            u_it++;
                            continue;
                        }
                        else
                        {
                            ns = new_ns;
                            update(pos_vertices, neg_vertices, deleted_vertices, cc, sim_graph);
                            Removed_vertices.push_back(vertexInfo);
                            for (auto &v : deleted_vertices)
                            {
                                cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
                            }

                            break;
                        }
                    }
                }
            }
            if (u_it == vec.end())
            {
                break;
            }
        }
    }

    if (ns >= thd1 && ns <= thd2)
    {
        Weight score = WEIGHT_MAX;
        for (int i = 0; i < cc.size(); i++)
        {
            int u = cc[i];

            for (int j = i + 1; j < cc.size(); j++)
            {
                int v = cc[j];

                score = min(score, sim_graph[u][v]);
            }
        }
        // cout << "score: " << score.numerator << " / " << score.denominator << endl;
        // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
        if (score > bottom_score)
        {
            bottom_score = score;
            // printGraphAsEdgeList(graph, cc);
            // if (bottom_score == Weight(1, 526))
            // // if (bottom_score == Weight(1, 192))
            // // if (bottom_score == Weight(1, 6336))
            // {
            //     // printGraphAsEdgeList(graph, cc);
            //     int a = 1;
            // }
            // cout << "bottom_score: " << bottom_score.numerator << "/" << bottom_score.denominator << endl;
            // if (bottom_score == Weight(1, 100))
            // {
            //     printGraphAsEdgeList(graph, cc);
            //     for (int i = 0; i < cc.size(); i++)
            //     {
            //         int u = cc[i];

            //         for (int j = i + 1; j < cc.size(); j++)
            //         {
            //             int v = cc[j];

            //             cout << u << "," << v << " " << sim_graph[u][v].numerator << "/" << sim_graph[u][v].denominator << endl;
            //         }
            //     }
            // }
            // if (bottom_score == Weight(1, 100))
            // {
            //     printGraphAsEdgeList(graph, cc);
            // }
            // remove vertex with score less than bottom_score and maintain the truss
            vector<VertexInfo> changed_edges;
            vector<set<int>> del_vertices;
            greedy_remove(changed_edges, del_vertices, cc, graph, sim_graph, bottom_score, k);
            // find connected components
            vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
            // calculate all added vertex-pair score
            // set<int> added_vertices;
            for (auto &cc : connectedComponents)
            {
                if (cc.size() <= 1)
                    continue;
                // repeat the process
                // printGraphAsEdgeList(graph, cc);
                // cout << "delete edges" << endl;
                special_shrink(special, graph, cc, bottom_score, seq2att, k, thd1, thd2, sim_graph, startTime);
            }
            // for (auto &it : del_vertices)
            // {
            //     for (int v : it)
            //     {
            //         graph.emplace(v, unordered_map<int, int>());
            //     }
            // }

            while (!changed_edges.empty())
            {
                auto &it = changed_edges.back();

                recoverTruss(it, graph);
                changed_edges.pop_back();
            }
        }
    }

    // recovery
    // if (score_id == score_edges.size())
    //     break;
    while (!Removed_vertices.empty())
    {
        VertexInfo v = Removed_vertices.back();
        // RecoverEdges(v.vertex, v.neighbors, v.changed_edges, graph);
        recoverTruss(v, graph);
        Removed_vertices.pop_back();
    }
}
// void shrink(std::unordered_map<int, std::unordered_map<int, int>> &graph, vector<int> &cc, Weight &bottom_score,
//             map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
//             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//             unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     Weight ns = NSize(cc, seq2att);
//     if (cc.size() == k)
//     {
//         if (ns >= thd1 && ns <= thd2)
//         {
//             Weight score = WEIGHT_MAX;
//             for (int i = 0; i < cc.size(); i++)
//             {
//                 int u = cc[i];

//                 for (int j = i + 1; j < cc.size(); j++)
//                 {
//                     int v = cc[j];

//                     score = min(score, sim_graph[u][v]);
//                 }
//             }
//             bottom_score = score;
//             // printGraphAsEdgeList(graph, cc);
//         }
//         return;
//     }
//     if (thd2 < 1 && ns == Weight(1, 1))
//     {
//         return;
//     }
//     if (thd1 > 0 && ns == Weight(0, 1))
//     {
//         return;
//     }
//     vector<VertexInfo> Removed_vertices;
//     // unordered_set<int> pos_edged_vertices;
//     // unordered_set<int> neg_edged_vertices;
//     unordered_map<int, double> pos_edged_vertices;
//     unordered_map<int, double> neg_edged_vertices;
//     for (int v : cc)
//     {
//         // if(!graph.count(src)) continue;
//         // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//         double weight = calcul_weight(graph, cc, v, sim_graph);
//         // Weight weight = calcul_weight(graph, cc, v, sim_graph);
//         if (seq2att[v].first)
//         {
//             pos_edged_vertices.emplace(v, weight);
//         }
//         else
//         {
//             // neg_edged_vertices.emplace(v, sup);
//             neg_edged_vertices.emplace(v, weight);
//         }
//     }
//     // double ns = NSize(edge_truss, seq2att);
//     // double ns = 0.0;
//     // for (int u : cc)
//     // {
//     //     if (!seq2att[u].first)
//     //     {
//     //         ns++;
//     //     }
//     // }
//     // ns /= cc.size();
//     while (cc.size() > k && (ns < thd1 || ns > thd2))
//     {
//         // bool has_vertex = true;

//         if (ns < thd1)
//         {
//             // remove the pos vertices
//             vector<pair<int, double>> vec(pos_edged_vertices.begin(), pos_edged_vertices.end());
//             sort(vec.begin(), vec.end(), compareByValue);
//             auto u_it = vec.begin();
//             for (; u_it != vec.end();)
//             {
//                 int u = u_it->first;
//                 unordered_set<int> deleted_vertices;
//                 unordered_map<int, unordered_map<int, int>> deleted_edges;
//                 // vector<Edge> changed_edges;
//                 Deleted_vertex(u, graph, cc, k, deleted_vertices, deleted_edges);
//                 if (deleted_vertices.empty())
//                 {
//                     u_it++;
//                     continue;
//                 }
//                 int neg_num = 0;
//                 for (int v : deleted_vertices)
//                 {
//                     if (!seq2att[v].first)
//                     {
//                         neg_num++;
//                     }
//                 }
//                 Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - deleted_vertices.size());
//                 if (new_ns < ns)
//                 {
//                     u_it++;
//                     pos_edged_vertices.erase(u);
//                     continue;
//                 }
//                 else
//                 {
//                     ns = new_ns;
//                     // update pos/neg edged vertices
//                     for (int u : deleted_vertices)
//                     {
//                         if (pos_edged_vertices.count(u))
//                         {
//                             pos_edged_vertices.erase(u);
//                         }
//                         else if (neg_edged_vertices.count(u))
//                         {
//                             neg_edged_vertices.erase(u);
//                         }
//                     }
//                     // unordered_map<int, double> pos_neighbors2delete_weight;
//                     // unordered_map<int, int> pos_neighbors2delete_count;
//                     // unordered_map<int, double> neg_neighbors2delete_weight;
//                     // unordered_map<int, int> neg_neighbors2delete_count;
//                     // for (int v : deleted_vertices)
//                     // {
//                     //     for (auto &w_s : graph[v])
//                     //     {
//                     //         int w = w_s.first;
//                     //         if (pos_edged_vertices.count(w))
//                     //         {
//                     //             Weight a = sim_graph[v][w] + 1;
//                     //             // int neighbor_count = graph[w].size();
//                     //             // double b = pos_edged_vertices[w] * neighbor_count;
//                     //             if(!pos_neighbors2delete_weight.count(w)){
//                     //                 pos_neighbors2delete_weight.emplace(w, 0.0);
//                     //                 pos_neighbors2delete_count.emplace(w, 0);
//                     //             }
//                     //             pos_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                     //             pos_neighbors2delete_count[w] += 1;
//                     //         }
//                     //         else if(neg_edged_vertices.count(w)){
//                     //             Weight a = sim_graph[v][w] + 1;
//                     //             // neg_edged_vertices[w] -= log(a.numerator) - log(a.denominator);
//                     //             if(!neg_neighbors2delete_weight.count(w)){
//                     //                 neg_neighbors2delete_weight.emplace(w, 0.0);
//                     //                 neg_neighbors2delete_count.emplace(w, 0);
//                     //             }
//                     //             neg_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                     //             neg_neighbors2delete_count[w] += 1;
//                     //         }
//                     //     }
//                     // }
//                     // for(auto &v_w : pos_neighbors2delete_weight){
//                     //     int v = v_w.first;
//                     //     double ori_weight = pos_edged_vertices[v] * graph[v].size();
//                     //     pos_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - pos_neighbors2delete_count[v]);
//                     // }
//                     // for(auto &v_w : neg_neighbors2delete_weight){
//                     //     int v = v_w.first;
//                     //     double ori_weight = neg_edged_vertices[v] * graph[v].size();
//                     //     neg_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - neg_neighbors2delete_count[v]);

//                     // }
//                     for (auto &v_s : pos_edged_vertices)
//                     {
//                         int v = v_s.first;
//                         double ori_w = pos_edged_vertices[v] * cc.size();
//                         for (int w : deleted_vertices)
//                         {
//                             Weight a = sim_graph[w][v] + 1;
//                             ori_w -= log(a.numerator) - log(a.denominator);
//                         }
//                         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                     }
//                     for (auto &v_s : neg_edged_vertices)
//                     {
//                         int v = v_s.first;
//                         double ori_w = neg_edged_vertices[v] * cc.size();
//                         for (int w : deleted_vertices)
//                         {
//                             Weight a = sim_graph[w][v] + 1;
//                             ori_w -= log(a.numerator) - log(a.denominator);
//                         }
//                         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                     }
//                     // for (auto &v_s : pos_edged_vertices)
//                     // {
//                     //     int v = v_s.first;
//                     //     for (int w : deleted_vertices)
//                     //     {
//                     //         Weight a = sim_graph[w][v] + 1;
//                     //         v_s.second -= log(a.numerator) - log(a.denominator);
//                     //     }
//                     // }
//                     // for (auto &v_s : neg_edged_vertices)
//                     // {
//                     //     int v = v_s.first;
//                     //     for (int w : deleted_vertices)
//                     //     {
//                     //         Weight a = sim_graph[w][v] + 1;
//                     //         v_s.second -= log(a.numerator) - log(a.denominator);
//                     //     }
//                     // }

//                     unordered_set<int> new_vertex;
//                     for (int u : deleted_vertices)
//                     {
//                         for (auto &v_sup : graph[u])
//                         {
//                             int v = v_sup.first;
//                             // graph[v].erase(u); // remove u from v's neighbors

//                             // check whether has new edged vertex
//                             if (!new_vertex.count(v) && !pos_edged_vertices.count(v) && !neg_edged_vertices.count(v))
//                             {
//                                 new_vertex.insert(v);
//                                 // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                             }
//                         }
//                     }
//                     // update graph
//                     for (int u : deleted_vertices)
//                     {
//                         unordered_map<int, int> neighbors = graph[u];
//                         vector<Edge> changed_edges;
//                         for (auto &v_sup : neighbors)
//                         {
//                             int v = v_sup.first;
//                             graph[v].erase(u);
//                             graph[u].erase(v);
//                             for (auto &it1 : graph[v])
//                             {
//                                 int w = it1.first;
//                                 if (graph[u].count(w))
//                                 {
//                                     --graph[v][w];
//                                     --graph[w][v];
//                                     changed_edges.push_back(pair<int, int>(v, w));
//                                 }
//                             }
//                         }

//                         cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
//                         graph.erase(u);
//                         VertexInfo removeV;
//                         removeV.vertex = u;
//                         removeV.neighbors = move(neighbors);
//                         removeV.changed_edges = move(changed_edges);
//                         Removed_vertices.push_back(removeV);
//                     }
//                     // printGraphAsEdgeList(graph, cc);
//                     // for (auto &u_v : deleted_edges)
//                     // {
//                     //     int u = u_v.first;
//                     //     for (auto &v_sup : u_v.second)
//                     //     {
//                     //         int v = v_sup.first;
//                     //         graph[u].erase(v);
//                     //     }
//                     // }

//                     for (auto &v : new_vertex)
//                     {
//                         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                     }

//                     // for (auto &u_v : changed_edges)
//                     // {
//                     //     int u = u_v.first;
//                     //     int v = u_v.second;
//                     //     graph[u][v]--;
//                     //     graph[v][u]--;
//                     // }

//                     break;
//                 }
//             }
//             if (u_it == vec.end())
//             {
//                 break;
//             }
//         }
//         else if (ns > thd2)
//         {
//             unordered_map<int, double> new_neg_edged_vertices;
//             for (auto &u_w : neg_edged_vertices)
//             {
//                 int u = u_w.first;
//                 int pos_neighbor_size = 0;
//                 for (auto &u_v : graph[u])
//                 {
//                     int v = u_v.first;
//                     if (seq2att[v].first)
//                     {
//                         pos_neighbor_size++;
//                     }
//                 }
//                 new_neg_edged_vertices.emplace(u, u_w.second + pos_neighbor_size);
//             }
//             // vector<pair<int, double>> vec(neg_edged_vertices.begin(), neg_edged_vertices.end());
//             vector<pair<int, double>> vec(new_neg_edged_vertices.begin(), new_neg_edged_vertices.end());
//             sort(vec.begin(), vec.end(), compareByValue);
//             VertexInfo vertexInfo;
//             auto u_it = vec.begin();
//             for (; u_it != vec.end();)
//             {
//                 int u = u_it->first;
//                 set<int> del_vertex;

//                 vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, del_vertex);

//                 del_vertex.insert(u);
//                 bool connected = false;
//                 for (int v : cc)
//                 {
//                     if (!del_vertex.count(v))
//                     {
//                         unordered_set<int> visited;
//                         DFS(graph, v, visited);

//                         if (visited.size() == cc.size() - del_vertex.size())
//                         {
//                             connected = true;
//                         }
//                         break;
//                     }
//                 }

//                 if (!connected)
//                 {
//                     for (int v : del_vertex)
//                     {
//                         if (v != u)
//                         {
//                             graph.emplace(v, unordered_map<int, int>());
//                         }
//                     }
//                     recoverTruss(vertexInfo, graph);
//                     u_it++;
//                     continue;
//                 }
//                 else
//                 {
//                     // for (auto &v : del_vertex)
//                     // {
//                     //     cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
//                     // }
//                     int neg_num = 0;
//                     for (int v : del_vertex)
//                     {
//                         if (!seq2att[v].first)
//                         {
//                             neg_num++;
//                         }
//                     }
//                     Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - del_vertex.size());
//                     if (new_ns > ns)
//                     {
//                         for (int v : del_vertex)
//                         {
//                             if (v != u)
//                             {
//                                 graph.emplace(v, unordered_map<int, int>());
//                             }
//                         }
//                         recoverTruss(vertexInfo, graph);
//                         u_it++;
//                         continue;
//                     }
//                     else
//                     {
//                         ns = new_ns;
//                         Removed_vertices.push_back(vertexInfo);
//                         for (auto &v : del_vertex)
//                         {
//                             cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
//                             if (pos_edged_vertices.count(v))
//                             {
//                                 pos_edged_vertices.erase(v);
//                             }
//                             else if (neg_edged_vertices.count(v))
//                             {
//                                 neg_edged_vertices.erase(v);
//                             }
//                         }
//                         printGraphAsEdgeList(graph, cc);
//                         // update pos/neg edged vertices

//                         unordered_set<int> new_vertex;
//                         for (auto &e : vertexInfo.sup_less_edges)
//                         {
//                             int u = e.first.first;
//                             int v = e.first.second;
//                             if (!new_vertex.count(v) && !del_vertex.count(v))
//                             {
//                                 new_vertex.insert(v);
//                             }
//                             if (!new_vertex.count(u) && !del_vertex.count(u))
//                             {
//                                 new_vertex.insert(u);
//                             }
//                         }

//                         for (auto &v : new_vertex)
//                         {
//                             if(pos_edged_vertices.count(v)){
//                                 double ori_w = pos_edged_vertices[v] * cc.size();
//                                 for (int w : del_vertex)
//                                 {
//                                     Weight a = sim_graph[w][v] + 1;
//                                     // ori_w -= log(a.numerator) - log(a.denominator);
//                                     ori_w -= a.numerator / a.denominator;
//                                 }
//                                 // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
//                                 pos_edged_vertices[v] = ori_w / (cc.size() - del_vertex.size());
//                             }
//                             else if(neg_edged_vertices.count(v)){
//                                 double ori_w = neg_edged_vertices[v] * cc.size();
//                                 for (int w : del_vertex)
//                                 {
//                                     Weight a = sim_graph[w][v] + 1;
//                                     // ori_w -= log(a.numerator) - log(a.denominator);
//                                     ori_w -= a.numerator / a.denominator;
//                                 }
//                                 // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
//                                 neg_edged_vertices[v] = ori_w / (cc.size() - del_vertex.size());
//                             }
//                         }

//                         break;
//                     }
//                 }

//                 // //////////////////////////////////
//                 // unordered_set<int> deleted_vertices;
//                 // // unordered_map<int, unordered_map<int, int>> deleted_edges;
//                 // unordered_map<int, unordered_map<int, int>> deleted_edges;
//                 // Deleted_vertex(u, graph, cc, k, deleted_vertices, deleted_edges);
//                 // if (deleted_vertices.empty())
//                 // {
//                 //     u_it++;
//                 //     continue;
//                 // }
//                 // int neg_num = 0;
//                 // for (int v : deleted_vertices)
//                 // {
//                 //     if (!seq2att[v].first)
//                 //     {
//                 //         neg_num++;
//                 //     }
//                 // }
//                 // Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - deleted_vertices.size());
//                 // if (new_ns > ns)
//                 // {
//                 //     u_it++;
//                 //     neg_edged_vertices.erase(u);
//                 //     continue;
//                 // }
//                 // else
//                 // {
//                 //     ns = new_ns;
//                 //     // update pos/neg edged vertices
//                 //     for (int u : deleted_vertices)
//                 //     {
//                 //         if (pos_edged_vertices.count(u))
//                 //         {
//                 //             pos_edged_vertices.erase(u);
//                 //         }
//                 //         else if (neg_edged_vertices.count(u))
//                 //         {
//                 //             neg_edged_vertices.erase(u);
//                 //         }
//                 //     }
//                 //     // unordered_map<int, double> pos_neighbors2delete_weight;
//                 //     // unordered_map<int, int> pos_neighbors2delete_count;
//                 //     // unordered_map<int, double> neg_neighbors2delete_weight;
//                 //     // unordered_map<int, int> neg_neighbors2delete_count;
//                 //     // for (int v : deleted_vertices)
//                 //     // {
//                 //     //     for (auto &w_s : graph[v])
//                 //     //     {
//                 //     //         int w = w_s.first;
//                 //     //         if (pos_edged_vertices.count(w))
//                 //     //         {
//                 //     //             Weight a = sim_graph[v][w] + 1;
//                 //     //             // int neighbor_count = graph[w].size();
//                 //     //             // double b = pos_edged_vertices[w] * neighbor_count;
//                 //     //             if(!pos_neighbors2delete_weight.count(w)){
//                 //     //                 pos_neighbors2delete_weight.emplace(w, 0.0);
//                 //     //                 pos_neighbors2delete_count.emplace(w, 0);
//                 //     //             }
//                 //     //             pos_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                 //     //             pos_neighbors2delete_count[w] += 1;
//                 //     //         }
//                 //     //         else if(neg_edged_vertices.count(w)){
//                 //     //             Weight a = sim_graph[v][w] + 1;
//                 //     //             // neg_edged_vertices[w] -= log(a.numerator) - log(a.denominator);
//                 //     //             if(!neg_neighbors2delete_weight.count(w)){
//                 //     //                 neg_neighbors2delete_weight.emplace(w, 0.0);
//                 //     //                 neg_neighbors2delete_count.emplace(w, 0);
//                 //     //             }
//                 //     //             neg_neighbors2delete_weight[w] += log(a.numerator) - log(a.denominator);
//                 //     //             neg_neighbors2delete_count[w] += 1;
//                 //     //         }
//                 //     //     }
//                 //     // }
//                 //     // for(auto &v_w : pos_neighbors2delete_weight){
//                 //     //     int v = v_w.first;
//                 //     //     double ori_weight = pos_edged_vertices[v] * graph[v].size();
//                 //     //     pos_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - pos_neighbors2delete_count[v]);
//                 //     // }
//                 //     // for(auto &v_w : neg_neighbors2delete_weight){
//                 //     //     int v = v_w.first;
//                 //     //     double ori_weight = neg_edged_vertices[v] * graph[v].size();
//                 //     //     neg_edged_vertices[v] = (ori_weight - v_w.second) / (graph[v].size() - neg_neighbors2delete_count[v]);

//                 //     // }

//                 //     for (auto &v_s : pos_edged_vertices)
//                 //     {
//                 //         int v = v_s.first;
//                 //         double ori_w = pos_edged_vertices[v] * cc.size();
//                 //         for (int w : deleted_vertices)
//                 //         {
//                 //             Weight a = sim_graph[w][v] + 1;
//                 //             // ori_w -= log(a.numerator) - log(a.denominator);
//                 //             ori_w -= a.numerator / a.denominator;
//                 //         }
//                 //         // v_s.second = log(ori_w) - log(cc.size() - deleted_vertices.size());
//                 //         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                 //     }
//                 //     for (auto &v_s : neg_edged_vertices)
//                 //     {
//                 //         int v = v_s.first;
//                 //         double ori_w = neg_edged_vertices[v] * cc.size();
//                 //         for (int w : deleted_vertices)
//                 //         {
//                 //             Weight a = sim_graph[w][v] + 1;
//                 //             ori_w -= a.numerator / a.denominator;
//                 //         }
//                 //         v_s.second = ori_w / (cc.size() - deleted_vertices.size());
//                 //     }

//                 //     unordered_set<int> new_vertex;
//                 //     for (int u : deleted_vertices)
//                 //     {
//                 //         for (auto &v_sup : graph[u])
//                 //         {
//                 //             int v = v_sup.first;
//                 //             // graph[v].erase(u); // remove u from v's neighbors

//                 //             // check whether has new edged vertex
//                 //             if (!deleted_vertices.count(v) && !new_vertex.count(v) && !pos_edged_vertices.count(v) && !neg_edged_vertices.count(v))
//                 //             {
//                 //                 new_vertex.insert(v);
//                 //                 // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                 //             }
//                 //         }
//                 //     }
//                 //     // update graph
//                 //     for (int u : deleted_vertices)
//                 //     {
//                 //         unordered_map<int, int> neighbors = graph[u];
//                 //         vector<Edge> changed_edges;
//                 //         for (auto &v_sup : neighbors)
//                 //         {
//                 //             int v = v_sup.first;
//                 //             graph[v].erase(u);
//                 //             graph[u].erase(v);
//                 //             for (auto &it1 : graph[v])
//                 //             {
//                 //                 int w = it1.first;
//                 //                 if (graph[u].count(w))
//                 //                 {
//                 //                     --graph[v][w];
//                 //                     --graph[w][v];
//                 //                     changed_edges.push_back(pair<int, int>(v, w));
//                 //                 }
//                 //             }
//                 //         }

//                 //         cc.erase(remove(cc.begin(), cc.end(), u), cc.end());
//                 //         graph.erase(u);
//                 //         VertexInfo removeV;
//                 //         removeV.vertex = u;
//                 //         removeV.neighbors = move(neighbors);
//                 //         removeV.changed_edges = move(changed_edges);
//                 //         Removed_vertices.push_back(removeV);
//                 //     }
//                 //     // printGraphAsEdgeList(graph, cc);
//                 //     // for (auto &u_v : deleted_edges)
//                 //     // {
//                 //     //     int u = u_v.first;
//                 //     //     for (auto &v_sup : u_v.second)
//                 //     //     {
//                 //     //         int v = v_sup.first;
//                 //     //         graph[u].erase(v);
//                 //     //     }
//                 //     // }

//                 //     for (auto &v : new_vertex)
//                 //     {
//                 //         isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                 //     }

//                 //     break;
//                 // }
//             }
//             if (u_it == vec.end())
//             {
//                 break;
//                 // // vector<int> targets = {3, 25, 71, 168};
//                 // // if (containsAll(cc, targets))
//                 // // {
//                 // // printGraphAsEdgeList(graph, cc);
//                 // // }
//                 // // Weight score = WEIGHT_MAX;
//                 // // for (int i = 0; i < cc.size(); i++)
//                 // // {
//                 // //     int u = cc[i];

//                 // //     for (int j = i + 1; j < cc.size(); j++)
//                 // //     {
//                 // //         int v = cc[j];
//                 // //         if(sim_graph[u][v] == Weight(1, 204)){
//                 // //             int a = 1;
//                 // //         }
//                 // //         score = min(score, sim_graph[u][v]);
//                 // //     }
//                 // // }
//                 // // break;
//                 // /*clear pos/neg set, find the non-cut vertex,
//                 // select the vertex don't destroy the opposite sides,
//                 // find the new pos/neg vertex from the deleted vertex*/
//                 // // set<int> cv;
//                 // // BCC(graph, cv);

//                 // // connect all pos vertex and expand to k-truss
//                 // // find the minimum degree neg vertex and remove it, check the remaining subgraph whether is connected
//                 // vector<pair<int, int>> neg_vertices;
//                 // for (auto &v : cc)
//                 // {
//                 //     if (!seq2att[v].first)
//                 //     {
//                 //         neg_vertices.push_back(make_pair(v, graph[v].size()));
//                 //     }
//                 // }
//                 // sort(neg_vertices.begin(), neg_vertices.end(),
//                 //      [](const pair<int, int> &a, const pair<int, int> &b)
//                 //      {
//                 //          return a.second < b.second;
//                 //      });
//                 // VertexInfo vertexInfo;
//                 // auto u_it = neg_vertices.begin();
//                 // // select a neg vertex
//                 // for (; u_it != neg_vertices.end();)
//                 // {
//                 //     int u = u_it->first;
//                 //     set<int> del_vertex;

//                 //     vertexInfo = MaintainTruss_DeleteVertex(u, graph, k, del_vertex);
//                 //     del_vertex.insert(u);
//                 //     bool connected = false;
//                 //     for (int v : cc)
//                 //     {
//                 //         if (!del_vertex.count(v))
//                 //         {
//                 //             unordered_set<int> visited;
//                 //             DFS(graph, v, visited);

//                 //             if (visited.size() == cc.size() - del_vertex.size())
//                 //             {
//                 //                 connected = true;
//                 //             }
//                 //             break;
//                 //         }
//                 //     }

//                 //     if (!connected)
//                 //     {
//                 //         for (int v : del_vertex)
//                 //         {
//                 //             if (v != u)
//                 //             {
//                 //                 graph.emplace(v, unordered_map<int, int>());
//                 //             }
//                 //         }
//                 //         recoverTruss(vertexInfo, graph);
//                 //         u_it++;
//                 //         continue;
//                 //     }
//                 //     else
//                 //     {
//                 //         // for (auto &v : del_vertex)
//                 //         // {
//                 //         //     cc.erase(remove(cc.begin(), cc.end(), v), cc.end());
//                 //         // }
//                 //         int neg_num = 0;
//                 //         for (int v : del_vertex)
//                 //         {
//                 //             if (!seq2att[v].first)
//                 //             {
//                 //                 neg_num++;
//                 //             }
//                 //         }
//                 //         Weight new_ns = (ns * cc.size() - neg_num) / (cc.size() - del_vertex.size());
//                 //         if (new_ns > ns)
//                 //         {
//                 //             for (int v : del_vertex)
//                 //             {
//                 //                 if (v != u)
//                 //                 {
//                 //                     graph.emplace(v, unordered_map<int, int>());
//                 //                 }
//                 //             }
//                 //             recoverTruss(vertexInfo, graph);
//                 //             u_it++;
//                 //             continue;
//                 //         }
//                 //         else
//                 //         {
//                 //             ns = new_ns;
//                 //             // update pos/neg edged vertices
//                 //             pos_edged_vertices.clear();
//                 //             neg_edged_vertices.clear();

//                 //             unordered_set<int> new_vertex;
//                 //             for (auto &e : vertexInfo.sup_less_edges)
//                 //             {
//                 //                 int u = e.first.first;
//                 //                 int v = e.first.second;
//                 //                 if (!new_vertex.count(v) && !del_vertex.count(v))
//                 //                 {
//                 //                     new_vertex.insert(v);
//                 //                 }
//                 //                 if (!new_vertex.count(u) && !del_vertex.count(u))
//                 //                 {
//                 //                     new_vertex.insert(u);
//                 //                 }
//                 //             }
//                 //             // for (int u : del_vertex)
//                 //             // {
//                 //             //     for (auto &v_sup : graph[u])
//                 //             //     {
//                 //             //         int v = v_sup.first;
//                 //             //         // graph[v].erase(u); // remove u from v's neighbors

//                 //             //         // check whether has new edged vertex
//                 //             //         if (!del_vertex.count(v) && !new_vertex.count(v))
//                 //             //         {
//                 //             //             new_vertex.insert(v);
//                 //             //             // isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                 //             //         }
//                 //             //     }
//                 //             // }
//                 //             for (auto &v : new_vertex)
//                 //             {
//                 //                 isEdgedVertex(graph, cc, v, seq2att, pos_edged_vertices, neg_edged_vertices, sim_graph);
//                 //             }

//                 //             break;
//                 //         }
//                 //     }
//                 // }
//                 // if (u_it == neg_vertices.end())
//                 // {
//                 //     break;
//                 // }
//             }
//         }
//         // // remove the neg vertices
//         // if (First_Neg_Edged_Vertex(graph, cc, k, ns,
//         //                            seq2att, Removed_vertices,
//         //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//         // {
//         //     break;
//         //     // if (First_Pos_Edged_Vertex(graph, cc, k, ns,
//         //     //                            seq2att, Removed_vertices,
//         //     //                            pos_edged_vertices, neg_edged_vertices, sim_graph))
//         //     // {
//         //     //     // has_vertex =  false;
//         //     //     break;
//         //     // }
//         // }
//     }

//     if (ns >= thd1 && ns <= thd2)
//     {
//         Weight score = WEIGHT_MAX;
//         for (int i = 0; i < cc.size(); i++)
//         {
//             int u = cc[i];

//             for (int j = i + 1; j < cc.size(); j++)
//             {
//                 int v = cc[j];

//                 score = min(score, sim_graph[u][v]);
//             }
//         }
//         // cout << "score: " << score.numerator << " / " << score.denominator << endl;
//         // cout << "ns: " << ns.numerator << " / " << ns.denominator << endl;
//         if (score > bottom_score)
//         {
//             bottom_score = score;
//             if (bottom_score == Weight(1, 526))
//             // if (bottom_score == Weight(1, 192))
//             // if (bottom_score == Weight(1, 6336))
//             {
//                 // printGraphAsEdgeList(graph, cc);
//                 int a = 1;
//             }
//             cout << "bottom_score: " << bottom_score.numerator << "/" << bottom_score.denominator << endl;
//             printGraphAsEdgeList(graph, cc);
//             // remove vertex with score less than bottom_score and maintain the truss
//             vector<VertexInfo> changed_edges;
//             vector<set<int>> del_vertices;
//             greedy_remove(changed_edges, del_vertices, cc, graph, sim_graph, bottom_score, k);
//             // find connected components
//             vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
//             // calculate all added vertex-pair score
//             // set<int> added_vertices;
//             for (auto &cc : connectedComponents)
//             {
//                 if (cc.size() <= 1)
//                     continue;
//                 // repeat the process
//                 // printGraphAsEdgeList(graph, cc);
//                 cout << "delete edges" << endl;
//                 shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//             }
//             for (auto &it : del_vertices)
//             {
//                 for (int v : it)
//                 {
//                     graph.emplace(v, unordered_map<int, int>());
//                 }
//             }

//             while (!changed_edges.empty())
//             {
//                 auto &it = changed_edges.back();

//                 recoverTruss(it, graph);
//                 changed_edges.pop_back();
//             }
//         }
//     }

//     // recovery
//     // if (score_id == score_edges.size())
//     //     break;
//     while (!Removed_vertices.empty())
//     {
//         VertexInfo v = Removed_vertices.back();
//         RecoverEdges(v.vertex, v.neighbors, v.changed_edges, graph);
//         Removed_vertices.pop_back();
//     }
// }
// Weight GNKCS(bool decom, map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
//              unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//              unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     Weight bottom_score = WEIGHT_ZERO;

//     // collect the maximum score edges

//     unordered_map<int, unordered_map<int, int>> graph; // k-truss
//     // unordered_map<int, unordered_map<int, int>> graph1;
//     // unordered_map<int, unordered_map<int, double>> sim_graph;
//     int score_id = 0;

//     for (auto &it_se : score_edges)
//     {
//         // if (it_se.first < bottom_score)
//         //     break;
//         score_id++;
//         // if (score_id == 11)
//         // if (score_id == 106)
//         if (score_id == 243)
//         {
//             int a = 1;
//         }
//         for (auto &it1 : it_se.second)
//         {
//             update_sup(graph, it1);
//         }
//         // cout << "add after" << endl;

//         // // find the connected k-truss, and the added vertices

//         vector<Edge> temp_sup_less_edges;
//         vector<Edge> sup_less_edges;
//         if (!decom)
//         {
//             for (auto &it : graph)
//             {
//                 int u = it.first;
//                 for (auto &it1 : it.second)
//                 {
//                     int v = it1.first;
//                     if (v > u)
//                         continue;
//                     int sup = it1.second;
//                     if (sup < k - 2)
//                     {
//                         temp_sup_less_edges.push_back(pair<int, int>(u, v));
//                     }
//                 }
//             }
//             MaintainKTruss_DeleteEdge(graph, sup_less_edges, temp_sup_less_edges, k);
//             for (auto it = graph.begin(); it != graph.end();)
//             {
//                 // 检查删除边后是否产生孤立点
//                 if (it->second.empty())
//                 {
//                     it = graph.erase(it);
//                 }
//                 else
//                 {
//                     it++;
//                 }
//             }
//         }
//         // cout << "ktruss" << endl;
//         // printGraphAsEdgeList(graph);
//         // GetKTrussInc(graph, sup_less_edges, it.second, k);
//         vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
//         // calculate all added vertex-pair score
//         // set<int> added_vertices;
//         for (auto &cc : connectedComponents)
//         {
//             if (cc.size() <= 1)
//                 continue;
//             // cout << "score id:" << score_id << endl;
//             // printGraphAsEdgeList(graph, cc);
//             for (int i = 0; i < cc.size(); i++)
//             {
//                 int u = cc[i];
//                 Weight u_score = seq2att[u].second.first;
//                 vector<int> u_att = seq2att[u].second.second;
//                 if (!sim_graph.count(u))
//                 {
//                     sim_graph.emplace(u, unordered_map<int, Weight>());
//                 }
//                 for (int j = i + 1; j < cc.size(); j++)
//                 {
//                     int v = cc[j];
//                     if (!sim_graph[u].count(v))
//                     {
//                         Weight v_score = seq2att[v].second.first;
//                         vector<int> v_att = seq2att[v].second.second;
//                         vector<int> com_att;
//                         set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//                         Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//                         Weight u_v = cal_score(u_score, v_score, u_v_sim);
//                         sim_graph[u].emplace(v, u_v);
//                         sim_graph[v].emplace(u, u_v);
//                     }
//                 }
//             }
//             shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);

//         }
//         if(bottom_score != WEIGHT_ZERO){
//             break;
//         }
//         // if (score_id == score_edges.size())
//         //     break;
//         for (auto &edge : sup_less_edges)
//         {
//             update_sup(graph, edge);
//         }
//     }
//     return bottom_score;
// }
bool containsAll(const std::vector<int> &vec, const std::vector<int> &targets)
{
    return std::all_of(targets.begin(), targets.end(), [&](int val)
                       { return std::find(vec.begin(), vec.end(), val) != vec.end(); });
}
// Weight GNKCS(bool decom, unordered_map<int, unordered_map<int, int>> graph,
//              unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//              unordered_map<int, unordered_map<int, Weight>> &sim_graph)
// {
//     Weight bottom_score = WEIGHT_ZERO;

//     vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
//     // 获得新k-truss和新边

//     // calculate all added vertex-pair score
//     // set<int> added_vertices;
//     int cc_id = 0;
//     for (auto &cc : connectedComponents)
//     {
//         if (cc.size() < k)
//             continue;

//         for (int i = 0; i < cc.size(); i++)
//         {
//             int u = cc[i];
//             Weight u_score = seq2att[u].second.first;
//             vector<int> u_att = seq2att[u].second.second;
//             if (!sim_graph.count(u))
//             {
//                 sim_graph.emplace(u, unordered_map<int, Weight>());
//             }
//             for (int j = i + 1; j < cc.size(); j++)
//             {
//                 int v = cc[j];
//                 if (!sim_graph[u].count(v))
//                 {
//                     if ((u == 6 && v == 266) || (u == 266 && v == 6))
//                     {
//                         int a = 1;
//                     }
//                     Weight v_score = seq2att[v].second.first;
//                     vector<int> v_att = seq2att[v].second.second;
//                     vector<int> com_att;
//                     set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//                     Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//                     Weight u_v = cal_score(u_score, v_score, u_v_sim);
//                     sim_graph[u].emplace(v, u_v);
//                     sim_graph[v].emplace(u, u_v);
//                 }
//             }
//         }
//         Weight score = WEIGHT_ZERO;

//         special_shrink(graph, cc, score, seq2att, k, thd1, thd2, sim_graph);
//         bottom_score = max(bottom_score, score);
//     }
//     return bottom_score;
// }
///////////////5-14/////////////////////////////
Weight GNKCS(bool special, bool decom, map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
             unordered_map<int, unordered_map<int, Weight>> &sim_graph,
             chrono::high_resolution_clock::time_point startTime)
{
    Weight bottom_score = WEIGHT_ZERO;

    // collect the maximum score edges

    unordered_map<int, unordered_map<int, int>> graph; // k-truss
    // unordered_map<int, unordered_map<int, int>> graph1;
    // unordered_map<int, unordered_map<int, double>> sim_graph;
    int score_id = 0;
    vector<vector<int>> cktruss;
    vector<Edge> sup_less_edges;
    int sn = score_edges.size();
    // for(auto &it : score_edges){
    //     sn += it.second.size();
    // }
    // cout << "sn: " << sn << endl;
    for (auto &it_se : score_edges)
    {
        if (it_se.first <= bottom_score)
            break;
        // score_id++;
        // cout << "score id:" << score_id << endl;
        // if (score_id == 11)
        // if (score_id == 106)

        // for (auto &it1 : it.second)
        // {
        //     int u = it1.first, v = it1.second;
        //     if (!graph1.count(u))
        //     {
        //         graph1.emplace(u, unordered_map<int, int>());
        //     }
        //     graph1[u].emplace(v, -1);
        //     if (!graph1.count(v))
        //     {
        //         graph1.emplace(v, unordered_map<int, int>());
        //     }
        //     graph1[v].emplace(u, -1);
        // }
        // for (auto &it1 : graph1)
        // {
        //     int v = it1.first;
        //     for (auto &it2 : it1.second)
        //     {
        //         int u = it2.first;
        //         if (v > u)
        //             continue;
        //         int sup = 0;
        //         for (auto &it3 : graph1[u])
        //         {
        //             int w = it3.first;
        //             if (graph1[v].count(w))
        //             {
        //                 sup++;
        //             }
        //         }
        //         graph1[u][v]=sup;
        //         graph1[v][u]=sup;
        //     }
        // }
        // cout << "add before" << endl;
        // printGraphAsEdgeList(graph);
        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
        for (auto &edge : sup_less_edges)
        {
            update_sup(graph, edge);
        }
        sup_less_edges.clear();
        // cout << "add after" << endl;

        // // find the connected k-truss, and the added vertices

        vector<Edge> temp_sup_less_edges;

        if (!decom)
        {
            for (auto &it : graph)
            {
                int u = it.first;
                for (auto &it1 : it.second)
                {
                    int v = it1.first;
                    if (v > u)
                        continue;
                    int sup = it1.second;
                    if (sup < k - 2)
                    {
                        temp_sup_less_edges.push_back(pair<int, int>(u, v));
                    }
                }
            }
            MaintainKTruss_DeleteEdge(graph, sup_less_edges, temp_sup_less_edges, k);
            for (auto it = graph.begin(); it != graph.end();)
            {
                // 检查删除边后是否产生孤立点
                if (it->second.empty())
                {
                    it = graph.erase(it);
                }
                else
                {
                    it++;
                }
            }
        }
        if (score_id == 1541)
        {
            // printGraphAsEdgeList(graph);
            int a = 1;
        }
        // cout << "ktruss" << endl;
        // printGraphAsEdgeList(graph);
        // GetKTrussInc(graph, sup_less_edges, it.second, k);
        vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
        // 获得新k-truss和新边

        // calculate all added vertex-pair score
        // set<int> added_vertices;
        int cc_id = 0;
        for (auto &cc : connectedComponents)
        {
            if (cc.size() < k)
                continue;
            bool old = false;
            cc_id++;
            // cout << "cc id:" << cc_id << endl;
            // if (score_id == 15 && cc_id == 1)
            // {
            //     printGraphAsEdgeList(graph, cc);
            // }
            for (auto &ktruss : cktruss)
            {
                if (ktruss.size() == cc.size())
                {
                    sort(cc.begin(), cc.end());
                    sort(ktruss.begin(), ktruss.end());
                    if (cc == ktruss)
                    {
                        old = true;
                        break;
                    }
                }
            }
            if (old)
                continue;
            // cout << "score id:" << score_id << endl;
            // for(int u : cc){
            //     cout << u << ", ";
            // }
            // cout << endl;
            // vector<int> targets1 = {314, 70, 87, 84, 306, 49, 17, 374, 375, 582, 18, 6, 5, 3, 304, 370, 266, 378, 189, 63, 193, 495, 720, 320, 364, 142, 153, 622, 123, 23, 24, 102, 59, 548, 12, 95, 369, 94};
            // sort(targets1.begin(), targets1.end());
            // if (containsAll(cc, targets1))
            // {
            //     // printGraphAsEdgeList(graph, cc);
            //     int a = 1;
            // }
            // vector<int> targets2 = {548, 378, 266, 23, 17, 49, 6};
            // sort(targets2.begin(), targets2.end());
            // if (containsAll(cc, targets2))
            // {
            //     printGraphAsEdgeList(graph, cc);
            //     int a = 1;
            // }
            // vector<int> targets = {3, 25, 71, 168};
            // // vector<int> targets = {3, 71, 168};
            // if (containsAll(cc, targets))
            // {
            //     printGraphAsEdgeList(graph, cc);
            //     int a = 1;
            // }
            // cout << "sim cal" << endl;
            for (int i = 0; i < cc.size(); i++)
            {
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = now - startTime;
                // cout << "Time used: " << duration.count() << "s" << endl;
                if (duration.count() > max_time)
                {
                    return bottom_score;
                }
                int u = cc[i];
                Weight u_score = seq2att[u].second.first;
                vector<int> u_att = seq2att[u].second.second;
                if (!sim_graph.count(u))
                {
                    sim_graph.emplace(u, unordered_map<int, Weight>());
                }
                for (int j = i + 1; j < cc.size(); j++)
                {
                    int v = cc[j];
                    if (!sim_graph[u].count(v))
                    {
                        if ((u == 6 && v == 266) || (u == 266 && v == 6))
                        {
                            int a = 1;
                        }
                        Weight v_score = seq2att[v].second.first;
                        vector<int> v_att = seq2att[v].second.second;
                        vector<int> com_att;
                        set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
                        Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
                        Weight u_v = cal_score(u_score, v_score, u_v_sim);
                        sim_graph[u].emplace(v, u_v);
                        sim_graph[v].emplace(u, u_v);
                    }
                }
            }
            Weight score = WEIGHT_ZERO;
            // shrink(graph, cc, score, score_edges, seq2att, k, thd1, thd2, sim_graph);
            // printGraphAsEdgeList(graph, cc);
            // for (int i : cc)
            // {
            //     cout << i << ',';
            // }
            // cout << endl;
            // vector<int> cc_copy = {168, 177, 3, 59, 23, 17, 378, 266, 102, 18, 5, 6, 36, 49, 513, 375, 25, 519, 71, 64, 101};
            // bool f = true;
            // for (int i : cc_copy)
            // {
            //     if (find(cc.begin(), cc.end(), i) == cc.end())
            //     {
            //         f = false;
            //         break;
            //     }
            // }
            // if (f)
            // {
            //     cout << "yes" << endl;
            // }
            // cout << "start shrink" << endl;
            special_shrink(special, graph, cc, score, seq2att, k, thd1, thd2, sim_graph, startTime);
            bottom_score = max(bottom_score, score);
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = now - startTime;
            // cout << "Time used: " << duration.count() << "s" << endl;
            if (duration.count() > max_time)
            {
                return bottom_score;
            }
        }

        // if (score_id == score_edges.size())
        //     break;
        // for (auto &edge : sup_less_edges)
        // {
        //     update_sup(graph, edge);
        // }
        cktruss.clear();
        cktruss.insert(cktruss.end(), connectedComponents.begin(), connectedComponents.end());
    }
    return bottom_score;

    // if meets the nsize constraint, then remove vertices, if the remaining graph not meets the nsize constraint, then add new edges.
    // until the score of added vertices is the minimum score, then stop
    // delete the edges whose min neighbor pair score is less than the maximum score, get the subgraph whose score is greater than the maximum score

    // shrink the subgraph to find a k-truss meeting the nsize constraint, the removed vertex can not consider the min nei score,
    // check whether exist the subgraph meeting the nsize constraint
    // based on the nsize to identify the sign of vertex, select the vertex that dont destroy the k-truss structure, then check the connection.
    // select the vertex that dont destroy the k-truss structure, then check the connection. edge vertices: (1) not edges be delete; (2) degree == trussness, such vertex must not be the cut vertex. Thus, iteratively remove edged vertex
}
// Weight GNKCS(bool decom, map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
//              unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
//              unordered_map<int, unordered_map<int, Weight>> &sim_graph, map<Edge, Weight> &edge2score)
// {
//     Weight bottom_score = WEIGHT_ZERO;

//     // collect the maximum score edges

//     unordered_map<int, unordered_map<int, int>> graph; // k-truss
//     // unordered_map<int, unordered_map<int, int>> graph1;
//     // unordered_map<int, unordered_map<int, double>> sim_graph;
//     int score_id = 0;
//     vector<vector<int>> cktruss;
//     vector<Edge> sup_less_edges;
//     for (auto &it_se : score_edges)
//     {
//         if (it_se.first < bottom_score)
//             break;
//         score_id++;
//         // cout << "score id:" << score_id << endl;
//         // if (score_id == 11)
//         // if (score_id == 106)
//         if (score_id == 243)
//         {
//             int a = 1;
//         }

//         for (auto &it1 : it_se.second)
//         {
//             update_sup(graph, it1);
//             // int u = it1.first;
//             // int v = it1.second;
//             // if (edge2score.count(it1))
//             // {
//             //     sim_graph[u][v] = edge2score[it1];
//             // }
//             // else
//             // {
//             //     sim_graph[u][v] = edge2score[Edge(v, u)];
//             // }
//         }
//         for (auto &edge : sup_less_edges)
//         {
//             update_sup(graph, edge);
//         }
//         sup_less_edges.clear();
//         // cout << "add after" << endl;

//         // // find the connected k-truss, and the added vertices

//         vector<Edge> temp_sup_less_edges;

//         if (!decom)
//         {
//             for (auto &it : graph)
//             {
//                 int u = it.first;
//                 for (auto &it1 : it.second)
//                 {
//                     int v = it1.first;
//                     if (v > u)
//                         continue;
//                     int sup = it1.second;
//                     if (sup < k - 2)
//                     {
//                         temp_sup_less_edges.push_back(pair<int, int>(u, v));
//                     }
//                 }
//             }
//             MaintainKTruss_DeleteEdge(graph, sup_less_edges, temp_sup_less_edges, k);
//             for (auto it = graph.begin(); it != graph.end();)
//             {
//                 // 检查删除边后是否产生孤立点
//                 if (it->second.empty())
//                 {
//                     it = graph.erase(it);
//                 }
//                 else
//                 {
//                     it++;
//                 }
//             }
//         }
//         // cout << "ktruss" << endl;
//         // printGraphAsEdgeList(graph);
//         // GetKTrussInc(graph, sup_less_edges, it.second, k);
//         vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
//         // 获得新k-truss和新边

//         // calculate all added vertex-pair score
//         // set<int> added_vertices;
//         int cc_id = 0;
//         for (auto &cc : connectedComponents)
//         {
//             if (cc.size() < k)
//                 continue;
//             bool old = false;
//             cc_id++;
//             // cout << "cc id:" << cc_id << endl;
//             // if (score_id == 15 && cc_id == 1)
//             // {
//             //     printGraphAsEdgeList(graph, cc);
//             // }
//             for (auto &ktruss : cktruss)
//             {
//                 if (ktruss.size() == cc.size())
//                 {
//                     sort(cc.begin(), cc.end());
//                     sort(ktruss.begin(), ktruss.end());
//                     if (cc == ktruss)
//                     {
//                         old = true;
//                         break;
//                     }
//                 }
//             }
//             if (old)
//                 continue;

//             for (int i = 0; i < cc.size(); i++)
//             {
//                 int u = cc[i];
//                 Weight u_score = seq2att[u].second.first;
//                 vector<int> u_att = seq2att[u].second.second;
//                 if (!sim_graph.count(u))
//                 {
//                     sim_graph.emplace(u, unordered_map<int, Weight>());
//                 }
//                 for (int j = i + 1; j < cc.size(); j++)
//                 {
//                     int v = cc[j];
//                     if (!sim_graph[u].count(v))
//                     {
//                         if ((u == 6 && v == 266) || (u == 266 && v == 6))
//                         {
//                             int a = 1;
//                         }
//                         Weight v_score = seq2att[v].second.first;
//                         vector<int> v_att = seq2att[v].second.second;
//                         vector<int> com_att;
//                         set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//                         Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//                         Weight u_v = cal_score(u_score, v_score, u_v_sim);
//                         sim_graph[u].emplace(v, u_v);
//                         sim_graph[v].emplace(u, u_v);
//                     }
//                 }
//             }
//             Weight score = WEIGHT_ZERO;
//             // shrink(graph, cc, score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//             special_shrink(graph, cc, score, score_edges, seq2att, k, thd1, thd2, sim_graph);
//             bottom_score = max(bottom_score, score);
//         }

//         cktruss.clear();
//         cktruss.insert(cktruss.end(), connectedComponents.begin(), connectedComponents.end());
//     }
//     return bottom_score;
// }
////////////////////////////////////////5-14/////////////////
Weight GNKCS(int &score_id, bool decom, map<Weight, vector<pair<int, int>>, CompareKeys> &score_edges,
             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, int k, Weight thd1, Weight thd2,
             unordered_map<int, unordered_map<int, int>> &graph,
             unordered_map<int, unordered_map<int, Weight>> &sim_graph)
{
    Weight bottom_score = WEIGHT_ZERO;

    // collect the maximum score edges

    // unordered_map<int, unordered_map<int, int>> graph; // k-truss
    // unordered_map<int, unordered_map<int, int>> graph1;
    // unordered_map<int, unordered_map<int, double>> sim_graph;
    // int score_id = 0;

    for (auto &it_se : score_edges)
    {

        if (bottom_score != WEIGHT_ZERO)
            break;
        score_id++;
        // score = it_se.first;
        // if (score_id == 11)
        // if (score_id == 106)
        if (score_id == 243)
        {
            int a = 1;
        }

        // printGraphAsEdgeList(graph);
        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
        // cout << "add after" << endl;

        // // find the connected k-truss, and the added vertices

        vector<Edge> temp_sup_less_edges;
        vector<Edge> sup_less_edges;
        if (!decom)
        {
            for (auto &it : graph)
            {
                int u = it.first;
                for (auto &it1 : it.second)
                {
                    int v = it1.first;
                    if (v > u)
                        continue;
                    int sup = it1.second;
                    if (sup < k - 2)
                    {
                        temp_sup_less_edges.push_back(pair<int, int>(u, v));
                    }
                }
            }
            MaintainKTruss_DeleteEdge(graph, sup_less_edges, temp_sup_less_edges, k);
            for (auto it = graph.begin(); it != graph.end();)
            {
                // 检查删除边后是否产生孤立点
                if (it->second.empty())
                {
                    it = graph.erase(it);
                }
                else
                {
                    it++;
                }
            }
        }
        // cout << "ktruss" << endl;
        // printGraphAsEdgeList(graph);
        // GetKTrussInc(graph, sup_less_edges, it.second, k);
        vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
        // calculate all added vertex-pair score
        // set<int> added_vertices;
        for (auto &cc : connectedComponents)
        {
            if (cc.size() <= 1)
                continue;
            // cout << "score id:" << score_id << endl;
            // printGraphAsEdgeList(graph, cc);
            for (int i = 0; i < cc.size(); i++)
            {
                int u = cc[i];
                Weight u_score = seq2att[u].second.first;
                vector<int> u_att = seq2att[u].second.second;
                if (!sim_graph.count(u))
                {
                    sim_graph.emplace(u, unordered_map<int, Weight>());
                }
                for (int j = i + 1; j < cc.size(); j++)
                {
                    int v = cc[j];
                    if (!sim_graph[u].count(v))
                    {
                        Weight v_score = seq2att[v].second.first;
                        vector<int> v_att = seq2att[v].second.second;
                        vector<int> com_att;
                        set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
                        Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
                        Weight u_v = cal_score(u_score, v_score, u_v_sim);
                        sim_graph[u].emplace(v, u_v);
                        sim_graph[v].emplace(u, u_v);
                    }
                }
            }
            shrink(graph, cc, bottom_score, score_edges, seq2att, k, thd1, thd2, sim_graph);
        }
        if (score_id == score_edges.size())
            break;
        for (auto &edge : sup_less_edges)
        {
            update_sup(graph, edge);
        }
    }
    return bottom_score;

    // if meets the nsize constraint, then remove vertices, if the remaining graph not meets the nsize constraint, then add new edges.
    // until the score of added vertices is the minimum score, then stop
    // delete the edges whose min neighbor pair score is less than the maximum score, get the subgraph whose score is greater than the maximum score

    // shrink the subgraph to find a k-truss meeting the nsize constraint, the removed vertex can not consider the min nei score,
    // check whether exist the subgraph meeting the nsize constraint
    // based on the nsize to identify the sign of vertex, select the vertex that dont destroy the k-truss structure, then check the connection.
    // select the vertex that dont destroy the k-truss structure, then check the connection. edge vertices: (1) not edges be delete; (2) degree == trussness, such vertex must not be the cut vertex. Thus, iteratively remove edged vertex
}

void method8(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
             vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
{
    vector<int> vert_sort;
    unordered_map<int, Weight> vert_score = VertexScore(seq2att, kct, C);
    vector<pair<int, Weight>> vec(vert_score.begin(), vert_score.end());

    // 按第二个元素（double 值）由小到大排序，因为C中元素是从后往前取的
    std::sort(vec.begin(), vec.end(), [](const std::pair<int, Weight> &a, const std::pair<int, Weight> &b)
              { return a.second < b.second; });
    for (auto it : vec)
    {
        vert_sort.push_back(it.first);
    }
    NaiveEnum_AttCntTrussNsizeSort(M, vert_sort, kct, seq2att, k, opt, sthd1, sthd2, result, vert_score);
}
void method9(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
             vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result,
             chrono::high_resolution_clock::time_point startTime)
{
    unordered_map<int, unordered_map<int, Weight>> graph_sim;
    for (int i = 0; i < C.size(); i++)
    {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - startTime;
        if (duration.count() > max_time)
        {
            // std::cout << "Timeout: Execution exceeded " << max_time << " second." << std::endl;
            break;
        }
        int u = C[i];
        Weight u_score = seq2att[u].second.first;
        vector<int> u_att = seq2att[u].second.second;
        for (int j = i + 1; j < C.size(); j++)
        {
            int v = C[j];
            // if (u < v)
            //     continue;
            Weight v_score = seq2att[v].second.first;
            vector<int> v_att = seq2att[v].second.second;
            vector<int> com_att;
            set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
            Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
            Weight u_v = cal_score(u_score, v_score, u_v_sim);
            graph_sim[u][v] = u_v;
            graph_sim[v][u] = u_v;
        }
    }
    // NaiveEnum_AttCntTrussNsizeSortNotConst(M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
    NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, kct, seq2att, k, opt, sthd1, sthd2, result, graph_sim, startTime);
}
// map<Weight, vector<pair<int, int>>, CompareKeys> calculate_score_edges(unordered_map<int, unordered_map<int, int>> &kct,
//                                                                        unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
//                                                                        vector<int> &C,
//                                                                        Weight &opt)
// {
//     map<Weight, vector<pair<int, int>>, CompareKeys> score_edges;
//     for (int u : C)
//     {
//         Weight u_score = seq2att[u].second.first;
//         vector<int> u_att = seq2att[u].second.second;
//         for (auto it2 : kct[u])
//         {
//             int v = it2.first;
//             if (u > v)
//                 continue;
//             Weight v_score = seq2att[v].second.first;
//             vector<int> v_att = seq2att[v].second.second;
//             vector<int> com_att;
//             set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//             Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//             Weight u_v = cal_score(u_score, v_score, u_v_sim);
//             if (u_v < opt || u_v == WEIGHT_ZERO)
//                 continue;
//             if (!score_edges.count(u_v))
//             {
//                 score_edges.emplace(u_v, vector<pair<int, int>>());
//                 // score_edges.emplace(u_v, unordered_map<int, set<int>>());
//             }
//             score_edges[u_v].push_back(pair<int, int>(u, v));
//         }
//     }
//     return move(score_edges);
// }
map<Weight, vector<pair<int, int>>, CompareKeys> calculate_score_edges(unordered_map<int, unordered_map<int, int>> &kct,
                                                                       unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    map<Weight, vector<pair<int, int>>, CompareKeys> score_edges;
    for (auto &it : kct)
    {
        int u = it.first;
        Weight u_score = seq2att[u].second.first;
        vector<int> u_att = seq2att[u].second.second;
        for (auto it2 : kct[u])
        {
            int v = it2.first;
            if (u > v)
                continue;
            Weight v_score = seq2att[v].second.first;
            vector<int> v_att = seq2att[v].second.second;
            vector<int> com_att;
            set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
            Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
            Weight u_v = cal_score(u_score, v_score, u_v_sim);
            // if (u_v == WEIGHT_ZERO)
            //     continue;
            if (!score_edges.count(u_v))
            {
                score_edges.emplace(u_v, vector<pair<int, int>>());
                // score_edges.emplace(u_v, unordered_map<int, set<int>>());
            }
            score_edges[u_v].push_back(pair<int, int>(u, v));
        }
    }
    return move(score_edges);
}
// void calculate_score_edges(unordered_map<int, unordered_map<int, int>> &kct,
//                            unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
//                            map<Edge, Weight> &edge_sim,
//                            map<Weight, vector<pair<int, int>>> &score_edges)
// {
//     for (auto &it : kct)
//     {
//         int u = it.first;
//         Weight u_score = seq2att[u].second.first;
//         vector<int> u_att = seq2att[u].second.second;
//         for (auto it2 : kct[u])
//         {
//             int v = it2.first;
//             if (u > v)
//                 continue;
//             Weight v_score = seq2att[v].second.first;
//             vector<int> v_att = seq2att[v].second.second;
//             vector<int> com_att;
//             set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//             Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//             Weight u_v = cal_score(u_score, v_score, u_v_sim);
//             // if (u_v == WEIGHT_ZERO)
//             //     continue;
//             if (!score_edges.count(u_v))
//             {
//                 score_edges.emplace(u_v, vector<pair<int, int>>());
//                 // score_edges.emplace(u_v, unordered_map<int, set<int>>());
//             }
//             Edge e(u, v);
//             score_edges[u_v].push_back(e);
//             edge_sim.emplace(e, u_v);
//         }
//     }
// }
// void method10(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
//               vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, vector<pair<Weight, vector<int>>> &result)
// {
//     auto start1 = std::chrono::high_resolution_clock::now();
//     // calculate all edges' score
//     map<Weight, vector<pair<int, int>>, CompareKeys> score_edges = calculate_score_edges(kct, seq2att, C, opt);
//     unordered_map<int, unordered_map<int, int>> graph;
//     unordered_map<int, unordered_map<int, Weight>> graph_sim;
//     Weight bottom_score = GNKCS(false, score_edges, seq2att, k, sthd1, sthd2, graph_sim);

//     cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
//     auto end1 = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> duration1 = end1 - start1;
//     cout << "greedy time: " << duration1.count() << "s" << endl;
//     // for (auto &it_se : score_edges)
//     // {
//     //     cout << "edge score: " << it_se.first.numerator << " / " << it_se.first.denominator << ", num = " << it_se.second.size() << endl;
//     // }
//     // delete the edge less than bottom_score
//     opt = max(opt, bottom_score);
//     // unordered_map<int, unordered_map<int, int>> graph;
//     for (auto &it_se : score_edges)
//     {
//         if (it_se.first < opt)
//         {
//             break;
//         }

//         for (auto &it1 : it_se.second)
//         {
//             update_sup(graph, it1);
//         }
//     }
//     // cout << "edge_num: " << edge_num << endl;
//     // // find the connected k-truss, and the added vertices
//     vector<Edge> temp_sup_less_edges;
//     for (auto &it : graph)
//     {
//         int u = it.first;
//         for (auto &it1 : it.second)
//         {
//             int v = it1.first;
//             if (v > u)
//                 continue;
//             int sup = it1.second;
//             if (sup < k - 2)
//             {
//                 temp_sup_less_edges.push_back(pair<int, int>(u, v));
//             }
//         }
//     }
//     MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
//     // printGraphAsEdgeList(graph);
//     for (auto it = graph.begin(); it != graph.end();)
//     {
//         // 检查删除边后是否产生孤立点
//         if (it->second.empty())
//         {
//             it = graph.erase(it);
//         }
//         else
//         {
//             it++;
//         }
//     }
//     vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
//     // for (auto &cc : connectedComponents)
//     // {
//     //     if (cc.size() <= 1)
//     //         continue;
//     //     // cout << "score id:" << score_id << endl;
//     //     // printGraphAsEdgeList(graph, cc);
//     //     for (int i = 0; i < cc.size(); i++)
//     //     {
//     //         int u = cc[i];
//     //         Weight u_score = seq2att[u].second.first;
//     //         vector<int> u_att = seq2att[u].second.second;
//     //         if (!graph_sim.count(u))
//     //         {
//     //             graph_sim.emplace(u, unordered_map<int, Weight>());
//     //         }
//     //         for (int j = i + 1; j < cc.size(); j++)
//     //         {
//     //             int v = cc[j];
//     //             if (!graph_sim[u].count(v))
//     //             {
//     //                 Weight v_score = seq2att[v].second.first;
//     //                 vector<int> v_att = seq2att[v].second.second;
//     //                 vector<int> com_att;
//     //                 set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//     //                 Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//     //                 Weight u_v = cal_score(u_score, v_score, u_v_sim);
//     //                 graph_sim[u].emplace(v, u_v);
//     //                 graph_sim[v].emplace(u, u_v);
//     //             }
//     //         }
//     //     }
//     // }

//     auto start = std::chrono::high_resolution_clock::now();
//     for (auto &cc : connectedComponents)
//     {
//         cout << "cc size: " << cc.size() << endl;
//         // NaiveEnum_AttCntTrussNsizeSortNotConst(M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
//         NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> duration = end - start;
//     cout << "search time: " << duration.count() << "s" << endl;
// }
void bronKerboschPivot(unordered_map<int, unordered_map<int, Weight>> &graph,
                       vector<int> &R,
                       vector<int> &P,
                       vector<int> &X,
                       vector<vector<int>> &cliques)
{
    if (P.empty() && X.empty())
    {
        cliques.push_back(R);
        return;
    }

    // 选择枢轴顶点 - 在P ∪ X中选择度最大的顶点
    unordered_set<int> PUX(P.begin(), P.end());
    PUX.insert(X.begin(), X.end());

    int pivot = *PUX.begin();
    size_t maxDegree = graph[pivot].size();
    for (int v : PUX)
    {
        if (graph[v].size() > maxDegree)
        {
            pivot = v;
            maxDegree = graph[v].size();
        }
    }

    // P \ N(pivot)
    unordered_set<int> pivotNeighbors;
    for (auto &v_w : graph[pivot])
    {
        pivotNeighbors.insert(v_w.first);
    }
    vector<int> P_diff_N;
    for (int v : P)
    {
        if (pivotNeighbors.find(v) == pivotNeighbors.end())
        {
            P_diff_N.push_back(v);
        }
    }

    for (int v : P_diff_N)
    {
        // 获取当前顶点v的邻居（关键修正点）
        unordered_set<int> neighbors;
        for (auto &u_w : graph[v])
        {
            neighbors.insert(u_w.first); // 这里应该是neighbors不是pivotNeighbors
        }

        vector<int> R_new = R;
        R_new.push_back(v);

        vector<int> P_new;
        vector<int> X_new;

        // P ∩ N(v)
        for (int u : P)
        {
            if (neighbors.find(u) != neighbors.end())
            {
                P_new.push_back(u);
            }
        }

        // X ∩ N(v)
        for (int u : X)
        {
            if (neighbors.find(u) != neighbors.end())
            {
                X_new.push_back(u);
            }
        }

        bronKerboschPivot(graph, R_new, P_new, X_new, cliques);

        // 将v从P移到X
        P.erase(remove(P.begin(), P.end(), v), P.end());
        X.push_back(v);
    }
}
void method101(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
               vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, vector<pair<Weight, vector<int>>> &result)
{
    // auto start1 = std::chrono::high_resolution_clock::now();
    // // calculate all edges' score
    // map<Weight, vector<pair<int, int>>, CompareKeys> score_edges = calculate_score_edges(kct, seq2att, C, opt);
    // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // int score_id = 0;
    // Weight bottom_score = GNKCS(score_id, false, score_edges, seq2att, k, sthd1, sthd2, graph, graph_sim);

    // cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;

    // // for (auto &it_se : score_edges)
    // // {
    // //     cout << "edge score: " << it_se.first.numerator << " / " << it_se.first.denominator << ", num = " << it_se.second.size() << endl;
    // // }
    // // delete the edge less than bottom_score
    // opt = max(opt, bottom_score);
    // // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_set<int> old_vertices;
    // for (auto &u_entry : graph_sim)
    // {
    //     int u = u_entry.first;
    //     old_vertices.insert(u);

    //     auto &neighbors = u_entry.second; // 避免重复查找
    //     for (auto it = neighbors.begin(); it != neighbors.end();)
    //     {
    //         int v = it->first;
    //         if (u < v)
    //         {
    //             ++it;
    //             continue;
    //         }

    //         Weight weight = it->second;

    //         if (weight < opt)
    //         {
    //             it = neighbors.erase(it); // 安全地删除并更新迭代器
    //             graph_sim[v].erase(u);
    //         }
    //         else
    //         {
    //             ++it;
    //         }
    //     }
    // }
    // for (auto &it_se : score_edges)
    // {
    //     if (it_se.first < opt)
    //     {
    //         break;
    //     }
    //     if (score_id-- != 0)
    //         continue;

    //     for (auto &it1 : it_se.second)
    //     {
    //         update_sup(graph, it1);
    //     }
    // }
    // // cout << "edge_num: " << edge_num << endl;
    // // // find the connected k-truss, and the added vertices
    // vector<Edge> temp_sup_less_edges;
    // for (auto &it : graph)
    // {
    //     int u = it.first;
    //     for (auto &it1 : it.second)
    //     {
    //         int v = it1.first;
    //         if (v > u)
    //             continue;
    //         int sup = it1.second;
    //         if (sup < k - 2)
    //         {
    //             temp_sup_less_edges.push_back(pair<int, int>(u, v));
    //         }
    //     }
    // }
    // MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // // printGraphAsEdgeList(graph);
    // for (auto it = graph.begin(); it != graph.end();)
    // {
    //     // 检查删除边后是否产生孤立点
    //     if (it->second.empty())
    //     {
    //         it = graph.erase(it);
    //     }
    //     else
    //     {
    //         it++;
    //     }
    // }
    // vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // int graphSize = 0;
    // for (auto &cc : connectedComponents)
    // {
    //     graphSize += cc.size();
    //     cout << "Connected Component size: " << cc.size() << endl;
    //     if (cc.size() <= 1)
    //         continue;
    //     // cout << "score id:" << score_id << endl;
    //     // printGraphAsEdgeList(graph, cc);
    //     vector<int> new_vertices;
    //     vector<int> old_cc_vertices;
    //     for (int u : cc)
    //     {
    //         if (old_vertices.count(u))
    //         {
    //             old_cc_vertices.push_back(u);
    //         }
    //         else
    //         {
    //             new_vertices.push_back(u);
    //         }
    //     }
    //     for (int i = 0; i < new_vertices.size(); i++)
    //     {
    //         int u = new_vertices[i];
    //         Weight u_score = seq2att[u].second.first;
    //         vector<int> u_att = seq2att[u].second.second;
    //         if (!graph_sim.count(u))
    //         {
    //             graph_sim.emplace(u, unordered_map<int, Weight>());
    //         }
    //         for (int j = i + 1; j < new_vertices.size(); j++)
    //         {
    //             int v = new_vertices[j];
    //             if (!graph_sim[u].count(v))
    //             {
    //                 Weight v_score = seq2att[v].second.first;
    //                 vector<int> v_att = seq2att[v].second.second;
    //                 vector<int> com_att;
    //                 set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
    //                 Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
    //                 Weight u_v = cal_score(u_score, v_score, u_v_sim);
    //                 if (u_v < opt)
    //                     continue;
    //                 graph_sim[u].emplace(v, u_v);
    //                 graph_sim[v].emplace(u, u_v);
    //             }
    //         }
    //         for (int j = 0; j < old_cc_vertices.size(); j++)
    //         {
    //             int v = old_cc_vertices[j];
    //             if (!graph_sim[u].count(v))
    //             {
    //                 Weight v_score = seq2att[v].second.first;
    //                 vector<int> v_att = seq2att[v].second.second;
    //                 vector<int> com_att;
    //                 set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
    //                 Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
    //                 Weight u_v = cal_score(u_score, v_score, u_v_sim);
    //                 if (u_v < opt)
    //                     continue;
    //                 graph_sim[u].emplace(v, u_v);
    //                 graph_sim[v].emplace(u, u_v);
    //             }
    //         }
    //     }
    // }
    // vector<vector<int>> cliques;
    // vector<int> R; // 当前团
    // vector<int> P; // 候选节点
    // vector<int> X; // 已处理节点

    // for (auto it = graph_sim.begin(); it != graph_sim.end();)
    // {
    //     if (it->second.size() == 0)
    //     {
    //         it = graph_sim.erase(it);
    //     }
    //     else
    //     {
    //         P.push_back(it->first);
    //         it++;
    //     }
    // }
    // // bronKerbosch(graph_sim, R, P, X, 1, 0, P[1].size(), 0);
    // bronKerboschPivot(graph_sim, R, P, X, cliques);
    // // unordered_set<int> cans;
    // for (auto &cc : cliques)
    // {
    //     if (cc.size() < k)
    //         continue;
    //     // cans.insert(cc.begin(), cc.end());
    //     // connected k-truss
    //     unordered_map<int, unordered_map<int, int>> sub_graph;
    //     for (int u : cc)
    //     {
    //         for (auto &it : graph[u])
    //         {
    //             int v = it.first;
    //             if (u < v && count(cc.begin(), cc.end(), v))
    //             {
    //                 Edge edge(u, v);
    //                 update_sup(sub_graph, edge);
    //             }
    //         }
    //     }
    //     vector<Edge> temp_sup_less_edges;
    //     for (auto &it : sub_graph)
    //     {
    //         int u = it.first;
    //         for (auto &it1 : it.second)
    //         {
    //             int v = it1.first;
    //             if (v > u)
    //                 continue;
    //             int sup = it1.second;
    //             if (sup < k - 2)
    //             {
    //                 temp_sup_less_edges.push_back(pair<int, int>(u, v));
    //             }
    //         }
    //     }
    //     MaintainKTruss_DeleteEdge(sub_graph, temp_sup_less_edges, k);
    //     // printGraphAsEdgeList(graph);
    //     for (auto it = sub_graph.begin(); it != sub_graph.end();)
    //     {
    //         // 检查删除边后是否产生孤立点
    //         if (it->second.empty())
    //         {
    //             it = sub_graph.erase(it);
    //         }
    //         else
    //         {
    //             it++;
    //         }
    //     }
    //     vector<vector<int>> connectedComponents = FindConnectedComponents(sub_graph);
    //     // greedy
    // }

    // auto end1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration1 = end1 - start1;
    // cout << "greedy time: " << duration1.count() << "s" << endl;
    // auto start = std::chrono::high_resolution_clock::now();
    // for (auto &cc : connectedComponents)
    // {
    //     cout << "cc size: " << cc.size() << endl;
    //     // NaiveEnum_AttCntTrussNsizeSortNotConst(M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
    //     NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
    // }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;
}
void method103(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
               vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, vector<pair<Weight, vector<int>>> &result)
{
    // auto start1 = std::chrono::high_resolution_clock::now();
    // int less_edges = 0;
    // // calculate all edges' score
    // map<Weight, vector<pair<int, int>>, CompareKeys> score_edges = calculate_score_edges(kct, seq2att, C, opt);
    // cout << "less_edges: " << less_edges << endl;
    // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // int score_id = 0;
    // Weight bottom_score = GNKCS(score_id, false, score_edges, seq2att, k, sthd1, sthd2, graph, graph_sim);

    // cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;

    // // for (auto &it_se : score_edges)
    // // {
    // //     cout << "edge score: " << it_se.first.numerator << " / " << it_se.first.denominator << ", num = " << it_se.second.size() << endl;
    // // }
    // // delete the edge less than bottom_score
    // opt = max(opt, bottom_score);
    // // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_set<int> old_vertices;
    // for (auto &u_entry : graph_sim)
    // {
    //     int u = u_entry.first;
    //     old_vertices.insert(u);
    // }
    // for (auto &it_se : score_edges)
    // {
    //     if (it_se.first < opt)
    //     {
    //         break;
    //     }
    //     if (score_id-- != 0)
    //         continue;

    //     for (auto &it1 : it_se.second)
    //     {
    //         update_sup(graph, it1);
    //     }
    // }
    // // cout << "edge_num: " << edge_num << endl;
    // // // find the connected k-truss, and the added vertices
    // vector<Edge> temp_sup_less_edges;
    // for (auto &it : graph)
    // {
    //     int u = it.first;
    //     for (auto &it1 : it.second)
    //     {
    //         int v = it1.first;
    //         if (v > u)
    //             continue;
    //         int sup = it1.second;
    //         if (sup < k - 2)
    //         {
    //             temp_sup_less_edges.push_back(pair<int, int>(u, v));
    //         }
    //     }
    // }
    // MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // // printGraphAsEdgeList(graph);
    // for (auto it = graph.begin(); it != graph.end();)
    // {
    //     // 检查删除边后是否产生孤立点
    //     if (it->second.empty())
    //     {
    //         it = graph.erase(it);
    //     }
    //     else
    //     {
    //         it++;
    //     }
    // }
    // vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // int graphSize = 0;
    // for (auto &cc : connectedComponents)
    // {
    //     graphSize += cc.size();
    //     cout << "Connected Component size: " << cc.size() << endl;
    //     if (cc.size() <= 1)
    //         continue;
    //     // cout << "score id:" << score_id << endl;
    //     // printGraphAsEdgeList(graph, cc);
    //     vector<int> new_vertices;
    //     vector<int> old_cc_vertices;
    //     for (int u : cc)
    //     {
    //         if (old_vertices.count(u))
    //         {
    //             old_cc_vertices.push_back(u);
    //         }
    //         else
    //         {
    //             new_vertices.push_back(u);
    //         }
    //     }
    //     for (int i = 0; i < new_vertices.size(); i++)
    //     {
    //         int u = new_vertices[i];
    //         Weight u_score = seq2att[u].second.first;
    //         vector<int> u_att = seq2att[u].second.second;
    //         if (!graph_sim.count(u))
    //         {
    //             graph_sim.emplace(u, unordered_map<int, Weight>());
    //         }
    //         for (int j = i + 1; j < new_vertices.size(); j++)
    //         {
    //             int v = new_vertices[j];
    //             if (!graph_sim[u].count(v))
    //             {
    //                 Weight v_score = seq2att[v].second.first;
    //                 vector<int> v_att = seq2att[v].second.second;
    //                 vector<int> com_att;
    //                 set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
    //                 Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
    //                 Weight u_v = cal_score(u_score, v_score, u_v_sim);
    //                 graph_sim[u].emplace(v, u_v);
    //                 graph_sim[v].emplace(u, u_v);
    //             }
    //         }
    //         for (int j = 0; j < old_cc_vertices.size(); j++)
    //         {
    //             int v = old_cc_vertices[j];
    //             if (!graph_sim[u].count(v))
    //             {
    //                 Weight v_score = seq2att[v].second.first;
    //                 vector<int> v_att = seq2att[v].second.second;
    //                 vector<int> com_att;
    //                 set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
    //                 Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
    //                 Weight u_v = cal_score(u_score, v_score, u_v_sim);
    //                 graph_sim[u].emplace(v, u_v);
    //                 graph_sim[v].emplace(u, u_v);
    //             }
    //         }
    //     }
    //     // remove vertex
    //     vector<VertexInfo> changed_edges;
    //     vector<set<int>> del_vertices;
    //     greedy_remove(changed_edges, del_vertices, cc, graph, graph_sim, bottom_score, k);
    //     // find connected components
    //     vector<vector<int>> connectedComponents = FindConnectedComponents(graph, cc);
    //     // calculate all added vertex-pair score
    //     // set<int> added_vertices;
    //     for (auto &cc : connectedComponents)
    //     {
    //         if (cc.size() <= 1)
    //             continue;
    //         // repeat the process
    //         // printGraphAsEdgeList(graph, cc);
    //         cout << "delete edges" << endl;
    //         shrink(graph, cc, bottom_score, score_edges, seq2att, k, sthd1, sthd2, graph_sim);
    //     }
    // }
    // cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;

    // auto end1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration1 = end1 - start1;
    // cout << "greedy time: " << duration1.count() << "s" << endl;
    // auto start = std::chrono::high_resolution_clock::now();
    // for (auto &cc : connectedComponents)
    // {
    //     cout << "cc size: " << cc.size() << endl;
    //     // NaiveEnum_AttCntTrussNsizeSortNotConst(M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
    //     NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
    // }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;
}
void method102(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
               vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, vector<pair<Weight, vector<int>>> &result)
{
    // auto start1 = std::chrono::high_resolution_clock::now();
    // int less_edges = 0;
    // // calculate all edges' score
    // map<Weight, vector<pair<int, int>>, CompareKeys> score_edges = calculate_score_edges(kct, seq2att, C, opt);
    // cout << "less_edges: " << less_edges << endl;
    // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // int score_id = 0;
    // Weight bottom_score = GNKCS(score_id, false, score_edges, seq2att, k, sthd1, sthd2, graph, graph_sim);

    // cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;

    // // for (auto &it_se : score_edges)
    // // {
    // //     cout << "edge score: " << it_se.first.numerator << " / " << it_se.first.denominator << ", num = " << it_se.second.size() << endl;
    // // }
    // // delete the edge less than bottom_score
    // opt = max(opt, bottom_score);
    // // unordered_map<int, unordered_map<int, int>> graph;
    // // unordered_set<int> old_vertices;
    // // for (auto &u_entry : graph_sim)
    // // {
    // //     int u = u_entry.first;
    // //     old_vertices.insert(u);

    // //     auto &neighbors = u_entry.second; // 避免重复查找
    // //     for (auto it = neighbors.begin(); it != neighbors.end();)
    // //     {
    // //         int v = it->first;
    // //         if (u < v)
    // //         {
    // //             ++it;
    // //             continue;
    // //         }

    // //         Weight weight = it->second;

    // //         if (weight < opt)
    // //         {
    // //             it = neighbors.erase(it); // 安全地删除并更新迭代器
    // //             graph_sim[v].erase(u);
    // //         }
    // //         else
    // //         {
    // //             ++it;
    // //         }
    // //     }
    // // }
    // for (auto &it_se : score_edges)
    // {
    //     if (it_se.first < opt)
    //     {
    //         break;
    //     }
    //     if (score_id-- != 0)
    //         continue;

    //     for (auto &it1 : it_se.second)
    //     {
    //         update_sup(graph, it1);
    //     }
    // }
    // // cout << "edge_num: " << edge_num << endl;
    // // // find the connected k-truss, and the added vertices
    // vector<Edge> temp_sup_less_edges;
    // for (auto &it : graph)
    // {
    //     int u = it.first;
    //     for (auto &it1 : it.second)
    //     {
    //         int v = it1.first;
    //         if (v > u)
    //             continue;
    //         int sup = it1.second;
    //         if (sup < k - 2)
    //         {
    //             temp_sup_less_edges.push_back(pair<int, int>(u, v));
    //         }
    //     }
    // }
    // MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // // printGraphAsEdgeList(graph);
    // for (auto it = graph.begin(); it != graph.end();)
    // {
    //     // 检查删除边后是否产生孤立点
    //     if (it->second.empty())
    //     {
    //         it = graph.erase(it);
    //     }
    //     else
    //     {
    //         it++;
    //     }
    // }
    // vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // for (auto &cc : connectedComponents)
    // {
    //     if (cc.size() <= 1)
    //         continue;
    //     // cout << "score id:" << score_id << endl;
    //     // printGraphAsEdgeList(graph, cc);
    //     for (int i = 0; i < cc.size(); i++)
    //     {
    //         int u = cc[i];
    //         Weight u_score = seq2att[u].second.first;
    //         vector<int> u_att = seq2att[u].second.second;
    //         if (!graph_sim.count(u))
    //         {
    //             graph_sim.emplace(u, unordered_map<int, Weight>());
    //         }
    //         for (int j = i + 1; j < cc.size(); j++)
    //         {
    //             int v = cc[j];
    //             if (!graph_sim[u].count(v))
    //             {
    //                 Weight v_score = seq2att[v].second.first;
    //                 vector<int> v_att = seq2att[v].second.second;
    //                 vector<int> com_att;
    //                 set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
    //                 Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
    //                 Weight u_v = cal_score(u_score, v_score, u_v_sim);
    //                 graph_sim[u].emplace(v, u_v);
    //                 graph_sim[v].emplace(u, u_v);
    //             }
    //         }
    //     }
    // }
    // auto end1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration1 = end1 - start1;
    // cout << "greedy time: " << duration1.count() << "s" << endl;
    // auto start = std::chrono::high_resolution_clock::now();
    // for (auto &cc : connectedComponents)
    // {
    //     cout << "cc size: " << cc.size() << endl;
    //     // NaiveEnum_AttCntTrussNsizeSortNotConst(M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
    //     NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
    // }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;
}
// void method12(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
//               vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, vector<pair<Weight, vector<int>>> &result)
// {
// auto start1 = std::chrono::high_resolution_clock::now();
// int less_edges = 0;
// // calculate all edges' score
// map<Weight, vector<pair<int, int>>, CompareKeys> score_edges;

// // map<double, unordered_map<int, set<int>>, CompareKeys> score_edges;
// for (int u : C)
// {
//     Weight u_score = seq2att[u].second.first;
//     vector<int> u_att = seq2att[u].second.second;
//     for (auto it2 : kct[u])
//     {
//         int v = it2.first;
//         if (u > v)
//             continue;
//         Weight v_score = seq2att[v].second.first;
//         vector<int> v_att = seq2att[v].second.second;
//         vector<int> com_att;
//         set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//         Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//         // if (u_v_sim < Weight(1, 3151))
//         // if (u_v_sim < Weight(1, 643))
//         // if (u_v_sim < Weight(1, 192))
//         if (u_v_sim < Weight(13, 3618))
//         {
//             less_edges++;
//         }
//         Weight u_v = cal_score(u_score, v_score, u_v_sim);
//         if (!score_edges.count(u_v))
//         {
//             score_edges.emplace(u_v, vector<pair<int, int>>());
//             // score_edges.emplace(u_v, unordered_map<int, set<int>>());
//         }
//         score_edges[u_v].push_back(pair<int, int>(u, v));
//     }
// }
// cout << "less_edges: " << less_edges << endl;
// unordered_map<int, unordered_map<int, Weight>> graph_sim;
// Weight bottom_score = GNKCS(false, score_edges, seq2att, k, sthd1, sthd2, graph, graph_sim);

// cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;

// // for (auto &it_se : score_edges)
// // {
// //     cout << "edge score: " << it_se.first.numerator << " / " << it_se.first.denominator << ", num = " << it_se.second.size() << endl;
// // }
// // delete the edge less than bottom_score
// opt = max(opt, bottom_score);
// unordered_map<int, unordered_map<int, int>> graph;
// int edge_num = 0;
// for (auto &it_se : score_edges)
// {
//     if (it_se.first < opt)
//     {
//         break;
//     }
//     for (auto &it1 : it_se.second)
//     {
//         edge_num++;
//         update_sup(graph, it1);
//     }
// }
// // cout << "edge_num: " << edge_num << endl;
// // // find the connected k-truss, and the added vertices
// vector<Edge> temp_sup_less_edges;
// for (auto &it : graph)
// {
//     int u = it.first;
//     for (auto &it1 : it.second)
//     {
//         int v = it1.first;
//         if (v > u)
//             continue;
//         int sup = it1.second;
//         if (sup < k - 2)
//         {
//             temp_sup_less_edges.push_back(pair<int, int>(u, v));
//         }
//     }
// }
// MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
// // printGraphAsEdgeList(graph);
// for (auto it = graph.begin(); it != graph.end();)
// {
//     // 检查删除边后是否产生孤立点
//     if (it->second.empty())
//     {
//         it = graph.erase(it);
//     }
//     else
//     {
//         it++;
//     }
// }
// vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
// auto end1 = std::chrono::high_resolution_clock::now();
// std::chrono::duration<double> duration1 = end1 - start1;
// cout << "greedy time: " << duration1.count() << "s" << endl;
// auto start = std::chrono::high_resolution_clock::now();
// for (auto &cc : connectedComponents)
// {
//     cout << "cc size: " << cc.size() << endl;
//     // NaiveEnum_AttCntTrussNsizeSortNotConst(M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
//     NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
// }
// auto end = std::chrono::high_resolution_clock::now();
// std::chrono::duration<double> duration = end - start;
// cout << "search time: " << duration.count() << "s" << endl;
// }
// void method11(unordered_map<int, unordered_map<int, int>> &kct, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
//               vector<int> &M, vector<int> &C, int k, Weight &opt, Weight &sthd1, Weight &sthd2, vector<pair<Weight, vector<int>>> &result)
// {
//     auto start1 = std::chrono::high_resolution_clock::now();
//     map<Weight, vector<Edge>> s_edge;
//     map<Edge, Weight> edge_sim;
//     unordered_map<int, unordered_map<int, int>> CC;
//     for (int u : C)
//     {
//         CC.emplace(u, kct[u]);
//     }
//     for (int u : C)
//     {
//         Weight u_score = seq2att[u].second.first;
//         vector<int> u_att = seq2att[u].second.second;
//         for (auto it2 : CC[u])
//         {
//             int v = it2.first;
//             if (u > v)
//                 continue;
//             Weight v_score = seq2att[v].second.first;
//             vector<int> v_att = seq2att[v].second.second;
//             vector<int> com_att;
//             set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
//             Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
//             Weight u_v = cal_score(u_score, v_score, u_v_sim);
//             if (!s_edge.count(u_v))
//             {
//                 s_edge.emplace(u_v, vector<Edge>());
//             }
//             Edge e(u, v);
//             s_edge[u_v].push_back(e);
//             edge_sim.emplace(e, u_v);
//         }
//     }
//     KScoreTrussDecomposition *ks = new KScoreTrussDecomposition(CC, s_edge, edge_sim, k);
//     CC.clear();
//     s_edge.clear();
//     edge_sim.clear();
//     map<Edge, Weight> truss = ks->edge2score;
//     // calculate all edges' score
//     map<Weight, vector<pair<int, int>>, CompareKeys> score_edges;
//     for (auto &it : truss)
//     {
//         Weight s = it.second;
//         if (!score_edges.count(s))
//         {
//             score_edges.emplace(s, vector<pair<int, int>>());
//         }
//         score_edges[s].push_back(it.first);
//     }
//     truss.clear();
//     // map<double, unordered_map<int, set<int>>, CompareKeys> score_edges;

//     unordered_map<int, unordered_map<int, Weight>> graph_sim;
//     Weight bottom_score = GNKCS(true, score_edges, seq2att, k, sthd1, sthd2, graph_sim);

//     cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;

//     // for (auto &it_se : score_edges)
//     // {
//     //     cout << "edge score: " << it_se.first.numerator << " / " << it_se.first.denominator << ", num = " << it_se.second.size() << endl;
//     // }
//     // delete the edge less than bottom_score
//     opt = max(opt, bottom_score);
//     unordered_map<int, unordered_map<int, int>> graph;
//     int edge_num = 0;
//     for (auto &it_se : score_edges)
//     {
//         if (it_se.first < opt)
//         {
//             break;
//         }
//         for (auto &it1 : it_se.second)
//         {
//             edge_num++;
//             update_sup(graph, it1);
//         }
//     }
//     // cout << "edge_num: " << edge_num << endl;
//     // // find the connected k-truss, and the added vertices
//     // vector<Edge> temp_sup_less_edges;
//     // for (auto &it : graph)
//     // {
//     //     int u = it.first;
//     //     for (auto &it1 : it.second)
//     //     {
//     //         int v = it1.first;
//     //         if (v > u)
//     //             continue;
//     //         int sup = it1.second;
//     //         if (sup < k - 2)
//     //         {
//     //             temp_sup_less_edges.push_back(pair<int, int>(u, v));
//     //         }
//     //     }
//     // }
//     // MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
//     // // printGraphAsEdgeList(graph);
//     // for (auto it = graph.begin(); it != graph.end();)
//     // {
//     //     // 检查删除边后是否产生孤立点
//     //     if (it->second.empty())
//     //     {
//     //         it = graph.erase(it);
//     //     }
//     //     else
//     //     {
//     //         it++;
//     //     }
//     // }
//     vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
//     // sort(connectedComponents.begin(), connectedComponents.end(),
//     //      [](const std::vector<int> &a, const std::vector<int> &b)
//     //      {
//     //          return a.size() < b.size();
//     //      });
//     auto end1 = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> duration1 = end1 - start1;
//     cout << "greedy time: " << duration1.count() << "s" << endl;
//     auto start = std::chrono::high_resolution_clock::now();
//     for (auto &cc : connectedComponents)
//     {
//         cout << "cc size: " << cc.size() << endl;
//         NaiveEnum_AttCntTrussNsizeSortNotConst(M, cc, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
//     }
//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> duration = end - start;
//     cout << "search time: " << duration.count() << "s" << endl;
// }
void printf_results(DataGraph &dg, map<Weight, vector<vector<int>>, greater<Weight>> &result, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    if(result.empty()) return;
    cout << "Weight: " << result.begin()->first.toString() << endl;
    auto &results = result.begin()->second;
    for (int i = 0; i < results.size(); i++)
    {
        for (int j = 0; j < results[i].size(); j++)
        {
            int v = results[i][j];
            int v_id = dg.seq2id[v];
            results[i][j] = v_id;
            // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
            cout << v_id << ": " << seq2att[v].first << " ";
        }
        cout << endl;
        // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    }
}
map<Weight, vector<vector<int>>, greater<Weight>> Enumerate(int methods, string dataset, vector<int> posKws, vector<int> negKws,
                                                            int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    Weight opt = WEIGHT_ZERO;
    map<Weight, vector<vector<int>>, greater<Weight>> result;
    // subgraph based on keywords
    // auto start_time = std::chrono::high_resolution_clock::now();
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(0, AttFileName, posKws, negKws);
    std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    // find connected ktruss
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    // find connected components

    vector<vector<int>> connectedComponents = FindConnectedComponents(kct);

    // vector<Weight> result_score;
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();
    int ccId = 0;
    for (auto &C : connectedComponents)
    {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - start_time;
        if (duration.count() > max_time)
        {
            // std::cout << "Timeout: Execution exceeded " << max_time << " second." << std::endl;
            break;
        }
        // cout << "Connected Component: " << ccId++ << endl;
        // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
        vector<int> M;
        // double M_score = numeric_limits<double>::max();
        auto start = std::chrono::high_resolution_clock::now();
        if (methods == 0)
            NaiveEnum(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        else if (methods == 1)
            NaiveEnum_Att(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        // NaiveEnum_Att(M, C, kct, seq2att, k, opt, sthd, result, M_score);
        else if (methods == 2)
            NaiveEnum_Truss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        else if (methods == 3)
            NaiveEnum_Cnt(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        else if (methods == 4)
            NaiveEnum_CntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        else if (methods == 5)
            NaiveEnum_AttTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        else if (methods == 6)
            NaiveEnum_AttCntTruss(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        else if (methods == 7)
            NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        else if (methods == 8)
        {
            method8(kct, seq2att, M, C, k, opt, sthd1, sthd2, result);
        }
        else if (methods == 9)
        {
            method9(kct, seq2att, M, C, k, opt, sthd1, sthd2, result, start_time);
        }
        else if (methods == 10)
        {

            // method10(kct, seq2att, M, C, k, opt, sthd1, sthd2, result);
            // method103(kct, seq2att, M, C, k, opt, sthd1, sthd2, result);

            // get k-truss
            // get connected components
            //  bottom_score = 0;

            // NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, bottom_score, sthd1, sthd2, result);
        }
        else if (methods == 11)
        {
            // method11(kct, seq2att, M, C, k, opt, sthd1, sthd2, result);
            // get k-truss
            // get connected components
            //  bottom_score = 0;

            // NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, bottom_score, sthd1, sthd2, result);
        }
        else if (methods == 12)
        {
            // method12(kct, seq2att, M, C, k, opt, sthd1, sthd2, result);
            // get k-truss
            // get connected components
            //  bottom_score = 0;

            // NaiveEnum_AttCntTrussNsize(M, C, kct, seq2att, k, bottom_score, sthd1, sthd2, result);
        }
        auto end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        // cout << "time: " << duration.count() << "s" << endl;
    }
    printf_results(dg, result, seq2att);
    // for (int i = 0; i < result.size(); i++)
    // {
    //     cout << "Weight: " << result[i].first.toString() << ", ";
    //     for (int v : result[i].second)
    //     {
    //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
    //         cout << v << ": " << seq2att[v].first << " ";
    //     }
    //     cout << endl;
    //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    // }
    return move(result);
}
map<Weight, vector<vector<int>>, greater<Weight>> NoSpecial(int methods, string dataset, vector<int> posKws, vector<int> negKws,
                                                            int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    Weight opt = WEIGHT_ZERO;
    map<Weight, vector<vector<int>>, greater<Weight>> result;
    // subgraph based on keywords
    // auto start = std::chrono::high_resolution_clock::now();
    // auto start_time = std::chrono::high_resolution_clock::now();
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(0, AttFileName, posKws, negKws);
    std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    // find connected ktruss
    // cout << "Truss Decomposition: " << endl;
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();

    // auto start1 = std::chrono::high_resolution_clock::now();
    // calculate all edges' score
    map<Weight, vector<pair<int, int>>, CompareKeys> score_edges = calculate_score_edges(kct, seq2att);
    unordered_map<int, unordered_map<int, int>> graph;
    unordered_map<int, unordered_map<int, Weight>> graph_sim;
    map<Weight, vector<pair<int, int>>, CompareKeys> edges;
    edges.emplace(Weight(0, 1), vector<pair<int, int>>());
    auto &eds = edges.begin()->second;
    for (auto &it : score_edges)
    {

        Weight score = it.first;
        for (auto &it1 : it.second)
        {
            eds.push_back(it1);
            int u = it1.first;
            int v = it1.second;
            graph_sim[u][v] = score;
            graph_sim[v][u] = score;
        }
    }
    // cout << "start greedy" << endl;
    Weight bottom_score = GNKCS(false, true, edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    // duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;

    // cout << "sim vertex: " << graph_sim.size() << endl;
    cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - start_time;
    cout << "greedy time: " << duration1.count() << "s" << endl;
    if (duration1.count() > max_time)
    {
        return result;
    }
    opt = max(opt, bottom_score);
    graph.clear();
    // unordered_map<int, unordered_map<int, int>> graph;
    for (auto &it_se : score_edges)
    {
        if (it_se.first < opt)
        {
            break;
        }

        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
    }
    // cout << "edge_num: " << edge_num << endl;
    // // find the connected k-truss, and the added vertices
    vector<Edge> temp_sup_less_edges;
    for (auto &it : graph)
    {
        int u = it.first;
        for (auto &it1 : it.second)
        {
            int v = it1.first;
            if (v > u)
                continue;
            int sup = it1.second;
            if (sup < k - 2)
            {
                temp_sup_less_edges.push_back(pair<int, int>(u, v));
            }
        }
    }
    MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // printGraphAsEdgeList(graph);
    for (auto it = graph.begin(); it != graph.end();)
    {
        // 检查删除边后是否产生孤立点
        if (it->second.empty())
        {
            it = graph.erase(it);
        }
        else
        {
            it++;
        }
    }
    vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // find connected components
    // vector<vector<int>> connectedComponents = FindConnectedComponents(kct);

    // vector<Weight> result_score;

    int ccId = 0;
    for (auto &C : connectedComponents)
    {
        // cout << "Connected Component: " << ccId++ << endl;
        // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
        vector<int> M;
        // double M_score = numeric_limits<double>::max();
        sort(C.begin(), C.end());
        NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
    }
    printf_results(dg, result, seq2att);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;

    // for (int i = 0; i < result.size(); i++)
    // {
    //     cout << "Weight: " << result[i].first.toString() << ", ";
    //     for (int v : result[i].second)
    //     {
    //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
    //         cout << v << ": " << seq2att[v].first << " ";
    //     }
    //     cout << endl;
    //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    // }
    delete td;
    return move(result);
}
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerate(int methods, string dataset, vector<int> posKws, vector<int> negKws,
                                                               int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    auto greedy_start_time = std::chrono::high_resolution_clock::now();
    Weight opt = WEIGHT_ZERO;
    map<Weight, vector<vector<int>>, greater<Weight>> result;
    // subgraph based on keywords
    // auto start = std::chrono::high_resolution_clock::now();
    // auto start_time = std::chrono::high_resolution_clock::now();
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(0, AttFileName, posKws, negKws);
    std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    // find connected ktruss
    // cout << "Truss Decomposition: " << endl;
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();

    // auto start1 = std::chrono::high_resolution_clock::now();
    // calculate all edges' score
    map<Weight, vector<pair<int, int>>, CompareKeys> score_edges = calculate_score_edges(kct, seq2att);
    unordered_map<int, unordered_map<int, int>> graph;
    unordered_map<int, unordered_map<int, Weight>> graph_sim;
    map<Weight, vector<pair<int, int>>, CompareKeys> edges;
    edges.emplace(Weight(1, 1), vector<pair<int, int>>());
    auto &eds = edges.begin()->second;
    for (auto &it : score_edges)
    {

        Weight score = it.first;
        for (auto &it1 : it.second)
        {
            eds.push_back(it1);
            int u = it1.first;
            int v = it1.second;
            graph_sim[u][v] = score;
            graph_sim[v][u] = score;
        }
    }
    // cout << "start greedy" << endl;
    Weight bottom_score = GNKCS(true, true, edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    // duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;

    // cout << "sim vertex: " << graph_sim.size() << endl;
    cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - greedy_start_time;
    cout << "greedy time: " << duration1.count() << "s" << endl;
    duration1 = end1 - start_time;
    if (duration1.count() > max_time)
    {
        return result;
    }
    opt = max(opt, bottom_score);
    graph.clear();
    // unordered_map<int, unordered_map<int, int>> graph;
    for (auto &it_se : score_edges)
    {
        if (it_se.first < opt)
        {
            break;
        }

        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
    }
    // cout << "edge_num: " << edge_num << endl;
    // // find the connected k-truss, and the added vertices
    vector<Edge> temp_sup_less_edges;
    for (auto &it : graph)
    {
        int u = it.first;
        for (auto &it1 : it.second)
        {
            int v = it1.first;
            if (v > u)
                continue;
            int sup = it1.second;
            if (sup < k - 2)
            {
                temp_sup_less_edges.push_back(pair<int, int>(u, v));
            }
        }
    }
    MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // printGraphAsEdgeList(graph);
    for (auto it = graph.begin(); it != graph.end();)
    {
        // 检查删除边后是否产生孤立点
        if (it->second.empty())
        {
            it = graph.erase(it);
        }
        else
        {
            it++;
        }
    }
    vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // find connected components
    // vector<vector<int>> connectedComponents = FindConnectedComponents(kct);

    // vector<Weight> result_score;

    int ccId = 0;
    for (auto &C : connectedComponents)
    {
        // cout << "Connected Component: " << ccId++ << endl;
        // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
        vector<int> M;
        // double M_score = numeric_limits<double>::max();
        sort(C.begin(), C.end());
        if (methods == 0)
            NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
        else if (methods == 1)
            NaiveEnum_AttTruss(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
        else if (methods == 2)
            NaiveEnum_AttTrussCnnt(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
        else if (methods == 3)
            NaiveEnum_AttTrussSize(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
    }
    printf_results(dg, result, seq2att);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;

    // for (int i = 0; i < result.size(); i++)
    // {
    //     cout << "Weight: " << result[i].first.toString() << ", ";
    //     for (int v : result[i].second)
    //     {
    //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
    //         cout << v << ": " << seq2att[v].first << " ";
    //     }
    //     cout << endl;
    //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    // }
    delete td;
    return move(result);
}
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerateIncG(int methods, string dataset, vector<int> posKws, vector<int> negKws,
                                                                   int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    auto greedy_start_time = std::chrono::high_resolution_clock::now();
    Weight opt = WEIGHT_ZERO;
    map<Weight, vector<vector<int>>, greater<Weight>> result;
    // subgraph based on keywords
    // auto start = std::chrono::high_resolution_clock::now();
    // auto start_time = std::chrono::high_resolution_clock::now();
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(0, AttFileName, posKws, negKws);
    std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    // find connected ktruss
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();

    // auto start1 = std::chrono::high_resolution_clock::now();
    // calculate all edges' score
    map<Weight, vector<pair<int, int>>, CompareKeys> score_edges = calculate_score_edges(kct, seq2att);
    // unordered_map<int, unordered_map<int, int>> graph;
    unordered_map<int, unordered_map<int, Weight>> graph_sim;
    for (auto &it : score_edges)
    {
        Weight score = it.first;
        for (auto &it1 : it.second)
        {
            int u = it1.first;
            int v = it1.second;
            graph_sim[u][v] = score;
            graph_sim[v][u] = score;
        }
    }

    Weight bottom_score = GNKCS(true, false, score_edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    // duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;

    // cout << "sim vertex: " << graph_sim.size() << endl;
    cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - greedy_start_time;
    cout << "greedy time: " << duration1.count() << "s" << endl;
    duration1 = end1 - start_time;
    if (duration1.count() > max_time)
    {
        return result;
    }
    opt = max(opt, bottom_score);
    unordered_map<int, unordered_map<int, int>> graph;
    for (auto &it_se : score_edges)
    {
        if (it_se.first < opt)
        {
            break;
        }

        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
    }
    // cout << "edge_num: " << edge_num << endl;
    // // find the connected k-truss, and the added vertices
    vector<Edge> temp_sup_less_edges;
    for (auto &it : graph)
    {
        int u = it.first;
        for (auto &it1 : it.second)
        {
            int v = it1.first;
            if (v > u)
                continue;
            int sup = it1.second;
            if (sup < k - 2)
            {
                temp_sup_less_edges.push_back(pair<int, int>(u, v));
            }
        }
    }
    MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // printGraphAsEdgeList(graph);
    for (auto it = graph.begin(); it != graph.end();)
    {
        // 检查删除边后是否产生孤立点
        if (it->second.empty())
        {
            it = graph.erase(it);
        }
        else
        {
            it++;
        }
    }
    vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // find connected components
    // vector<vector<int>> connectedComponents = FindConnectedComponents(kct);

    // vector<Weight> result_score;

    int ccId = 0;
    for (auto &C : connectedComponents)
    {
        // cout << "Connected Component: " << ccId++ << endl;
        // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
        vector<int> M;
        // double M_score = numeric_limits<double>::max();
        sort(C.begin(), C.end());
        NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
    }
    printf_results(dg, result, seq2att);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;

    // for (int i = 0; i < result.size(); i++)
    // {
    //     cout << "Weight: " << result[i].first.toString() << ", ";
    //     for (int v : result[i].second)
    //     {
    //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
    //         cout << v << ": " << seq2att[v].first << " ";
    //     }
    //     cout << endl;
    //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    // }
    delete td;
    return move(result);
}
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerate1(int methods, string dataset, vector<int> posKws, vector<int> negKws,
                                                                int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    auto greedy_start_time = std::chrono::high_resolution_clock::now();
    Weight opt = WEIGHT_ZERO;
    map<Weight, vector<vector<int>>, greater<Weight>> result;
    // subgraph based on keywords
    // auto start = std::chrono::high_resolution_clock::now();
    // auto start_time = std::chrono::high_resolution_clock::now();
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(0, AttFileName, posKws, negKws);
    std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    // find connected ktruss
    // cout << "TrussDecomposition" << endl;
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();

    // auto start1 = std::chrono::high_resolution_clock::now();
    // calculate all edges' score
    map<Edge, Weight> edge_sim;
    map<Weight, vector<pair<int, int>>> score_edges;
    calculate_score_edges(kct, seq2att, edge_sim, score_edges);
    unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // for (auto &it : score_edges)
    // {
    //     Weight score = it.first;
    //     for (auto &it1 : it.second)
    //     {
    //         int u = it1.first;
    //         int v = it1.second;
    //         graph_sim[u][v] = score;
    //         graph_sim[v][u] = score;
    //     }
    // }
    // cout << "KScoreTrussDecomposition" << endl;
    KScoreTrussDecomposition *ks = new KScoreTrussDecomposition(kct, score_edges, edge_sim, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    edge_sim.clear();
    score_edges.clear();
    map<Edge, Weight> truss = ks->edge2score;
    // calculate all edges' score
    map<Weight, vector<pair<int, int>>, CompareKeys> s_truss_edges;
    for (auto &it : truss)
    {
        Weight s = it.second;
        if (!s_truss_edges.count(s))
        {
            s_truss_edges.emplace(s, vector<pair<int, int>>());
        }
        s_truss_edges[s].push_back(it.first);
    }
    truss.clear();
    // map<double, unordered_map<int, set<int>>, CompareKeys> score_edges;

    // cout << "GNKCS start" << endl;
    // Weight bottom_score = GNKCS(methods != 0, true, s_truss_edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    Weight bottom_score = GNKCS(true, true, s_truss_edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    // duration = std::chrono::high_resolution_clock::now() - start_time;
    // // cout << "Time used: " << duration.count() << "s" << endl;
    // if (duration.count() > max_time)
    // {
    //     return result;
    // }
    // cout << "sim vertex: " << graph_sim.size() << endl;
    // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // Weight bottom_score = GNKCS(false, score_edges, seq2att, k, sthd1, sthd2, graph_sim);

    cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - greedy_start_time;
    cout << "greedy time: " << duration1.count() << "s" << endl;
    duration1 = end1 - start_time;
    if (duration1.count() > max_time)
    {
        return result;
    }
    opt = max(opt, bottom_score);
    unordered_map<int, unordered_map<int, int>> graph;
    for (auto &it_se : s_truss_edges)
    {
        if (it_se.first < opt)
        {
            break;
        }

        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
    }
    // cout << "edge_num: " << edge_num << endl;
    // // find the connected k-truss, and the added vertices
    // vector<Edge> temp_sup_less_edges;
    // for (auto &it : graph)
    // {
    //     int u = it.first;
    //     for (auto &it1 : it.second)
    //     {
    //         int v = it1.first;
    //         if (v > u)
    //             continue;
    //         int sup = it1.second;
    //         if (sup < k - 2)
    //         {
    //             temp_sup_less_edges.push_back(pair<int, int>(u, v));
    //         }
    //     }
    // }
    // MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // // printGraphAsEdgeList(graph);
    // for (auto it = graph.begin(); it != graph.end();)
    // {
    //     // 检查删除边后是否产生孤立点
    //     if (it->second.empty())
    //     {
    //         it = graph.erase(it);
    //     }
    //     else
    //     {
    //         it++;
    //     }
    // }
    vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // find connected components
    // vector<vector<int>> connectedComponents = FindConnectedComponents(kct);

    // vector<Weight> result_score;

    int ccId = 0;
    for (auto &C : connectedComponents)
    {
        // cout << "Connected Component: " << ccId++ << endl;
        // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
        vector<int> M;
        // double M_score = numeric_limits<double>::max();
        /*for compare with another method, there fixed the order of traversed  vertices.*/
        sort(C.begin(), C.end());
        NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
    }
    printf_results(dg, result, seq2att);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;

    // for (int i = 0; i < result.size(); i++)
    // {
    //     cout << "Weight: " << result[i].first.toString() << ", ";
    //     for (int v : result[i].second)
    //     {
    //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
    //         cout << v << ": " << seq2att[v].first << " ";
    //     }
    //     cout << endl;
    //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    // }
    delete td;
    delete ks;
    return move(result);
}
map<Weight, vector<vector<int>>, greater<Weight>> AvdEnumerate2(int methods, string dataset, vector<int> posKws, vector<int> negKws,
                                                                int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    auto greedy_start_time = std::chrono::high_resolution_clock::now();
    Weight opt = WEIGHT_ZERO;
    map<Weight, vector<vector<int>>, greater<Weight>> result;
    // subgraph based on keywords
    // auto start = std::chrono::high_resolution_clock::now();
    // auto start_time = std::chrono::high_resolution_clock::now();
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(0, AttFileName, posKws, negKws);
    std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    // find connected ktruss
    // cout << "TrussDecomposition" << endl;
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();

    // auto start1 = std::chrono::high_resolution_clock::now();
    // calculate all edges' score
    map<Edge, Weight> edge_sim;
    map<Weight, vector<pair<int, int>>> score_edges;
    calculate_score_edges(kct, seq2att, edge_sim, score_edges);
    unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // for (auto &it : score_edges)
    // {
    //     Weight score = it.first;
    //     for (auto &it1 : it.second)
    //     {
    //         int u = it1.first;
    //         int v = it1.second;
    //         graph_sim[u][v] = score;
    //         graph_sim[v][u] = score;
    //     }
    // }
    // cout << "KScoreTrussDecomposition" << endl;
    KScoreTrussDecomposition *ks = new KScoreTrussDecomposition(kct, score_edges, edge_sim, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    edge_sim.clear();
    score_edges.clear();
    map<Edge, Weight> truss = ks->edge2score;
    // calculate all edges' score
    map<Weight, vector<pair<int, int>>, CompareKeys> s_truss_edges;
    for (auto &it : truss)
    {
        Weight s = it.second;
        if (!s_truss_edges.count(s))
        {
            s_truss_edges.emplace(s, vector<pair<int, int>>());
        }
        s_truss_edges[s].push_back(it.first);
    }
    truss.clear();
    // map<double, unordered_map<int, set<int>>, CompareKeys> score_edges;

    // cout << "GNKCS start" << endl;
    // Weight bottom_score = GNKCS(methods != 0, true, s_truss_edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    Weight bottom_score = GNKCS(true, true, s_truss_edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    // duration = std::chrono::high_resolution_clock::now() - start_time;
    // // cout << "Time used: " << duration.count() << "s" << endl;
    // if (duration.count() > max_time)
    // {
    //     return result;
    // }
    // cout << "sim vertex: " << graph_sim.size() << endl;
    // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // Weight bottom_score = GNKCS(false, score_edges, seq2att, k, sthd1, sthd2, graph_sim);

    cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - greedy_start_time;
    cout << "greedy time: " << duration1.count() << "s" << endl;
    duration1 = end1 - start_time;
    if (duration1.count() > max_time)
    {
        return result;
    }
    opt = max(opt, bottom_score);
    unordered_map<int, unordered_map<int, int>> graph;
    for (auto &it_se : s_truss_edges)
    {
        if (it_se.first < opt)
        {
            break;
        }

        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
    }
    // cout << "edge_num: " << edge_num << endl;
    // // find the connected k-truss, and the added vertices
    // vector<Edge> temp_sup_less_edges;
    // for (auto &it : graph)
    // {
    //     int u = it.first;
    //     for (auto &it1 : it.second)
    //     {
    //         int v = it1.first;
    //         if (v > u)
    //             continue;
    //         int sup = it1.second;
    //         if (sup < k - 2)
    //         {
    //             temp_sup_less_edges.push_back(pair<int, int>(u, v));
    //         }
    //     }
    // }
    // MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // // printGraphAsEdgeList(graph);
    // for (auto it = graph.begin(); it != graph.end();)
    // {
    //     // 检查删除边后是否产生孤立点
    //     if (it->second.empty())
    //     {
    //         it = graph.erase(it);
    //     }
    //     else
    //     {
    //         it++;
    //     }
    // }
    vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // find connected components
    // vector<vector<int>> connectedComponents = FindConnectedComponents(kct);

    // vector<Weight> result_score;

    int ccId = 0;
    for (auto &C : connectedComponents)
    {
        // cout << "Connected Component: " << ccId++ << endl;
        // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
        vector<int> M;
        // double M_score = numeric_limits<double>::max();
        /*for compare with another method, there fixed the order of traversed  vertices.*/
        sort(C.begin(), C.end());
        NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
    }
    printf_results(dg, result, seq2att);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;

    // for (int i = 0; i < result.size(); i++)
    // {
    //     cout << "Weight: " << result[i].first.toString() << ", ";
    //     for (int v : result[i].second)
    //     {
    //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
    //         cout << v << ": " << seq2att[v].first << " ";
    //     }
    //     cout << endl;
    //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    // }
    delete td;
    delete ks;
    return move(result);
}
// vector<pair<Weight, vector<int>>> AvdEnumerate2(int methods, string dataset, vector<int> posKws, vector<int> negKws, int k, Weight sthd1, Weight sthd2, Weight &opt)
// {
// unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(dataset, posKws, negKws);
// DataGraph dg = LoadGraph(dataset, id2att);
// // find connected ktruss
// TrussDecomposition *td = new TrussDecomposition(dg, k);
// unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
// unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
// for (auto &it : kct)
// {
//     int idx = dg.seq2id[it.first];
//     seq2att[it.first] = id2att[idx];
// }
// id2att.clear();

// auto start1 = std::chrono::high_resolution_clock::now();
// // calculate all edges' score
// map<Edge, Weight> edge_sim;

// map<Weight, vector<pair<int, int>>> score_edges;
// calculate_score_edges(kct, seq2att, edge_sim, score_edges);
// KScoreTrussDecomposition *ks = new KScoreTrussDecomposition(kct, score_edges, edge_sim, k);

// edge_sim.clear();
// score_edges.clear();
// map<Edge, Weight> truss = ks->edge2score;
// // calculate all edges' score
// map<Weight, vector<pair<int, int>>, CompareKeys> s_truss_edges;
// for (auto &it : truss)
// {
//     Weight s = it.second;
//     if (!s_truss_edges.count(s))
//     {
//         s_truss_edges.emplace(s, vector<pair<int, int>>());
//     }
//     s_truss_edges[s].push_back(it.first);
// }
// truss.clear();
// // map<double, unordered_map<int, set<int>>, CompareKeys> score_edges;

// unordered_map<int, unordered_map<int, Weight>> graph_sim;
// // for (auto &it : s_truss_edges)
// // {
// //     Weight score = it.first;
// //     for (auto &it1 : it.second)
// //     {
// //         int u = it1.first;
// //         int v = it1.second;
// //         graph_sim[u][v] = score;
// //         graph_sim[v][u] = score;
// //     }
// // }
// Weight bottom_score = GNKCS(true, s_truss_edges, seq2att, k, sthd1, sthd2, graph_sim);
// cout << "sim vertex: " << graph_sim.size() << endl;
// // unordered_map<int, unordered_map<int, int>> graph;
// // unordered_map<int, unordered_map<int, Weight>> graph_sim;
// // Weight bottom_score = GNKCS(false, score_edges, seq2att, k, sthd1, sthd2, graph_sim);

// cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
// auto end1 = std::chrono::high_resolution_clock::now();
// std::chrono::duration<double> duration1 = end1 - start1;
// cout << "greedy time: " << duration1.count() << "s" << endl;
// opt = max(opt, bottom_score);
// unordered_map<int, unordered_map<int, int>> graph;
// for (auto &it_se : s_truss_edges)
// {
//     if (it_se.first < opt)
//     {
//         break;
//     }

//     for (auto &it1 : it_se.second)
//     {
//         update_sup(graph, it1);
//     }
// }
// // cout << "edge_num: " << edge_num << endl;
// // // find the connected k-truss, and the added vertices
// // vector<Edge> temp_sup_less_edges;
// // for (auto &it : graph)
// // {
// //     int u = it.first;
// //     for (auto &it1 : it.second)
// //     {
// //         int v = it1.first;
// //         if (v > u)
// //             continue;
// //         int sup = it1.second;
// //         if (sup < k - 2)
// //         {
// //             temp_sup_less_edges.push_back(pair<int, int>(u, v));
// //         }
// //     }
// // }
// // MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
// // // printGraphAsEdgeList(graph);
// // for (auto it = graph.begin(); it != graph.end();)
// // {
// //     // 检查删除边后是否产生孤立点
// //     if (it->second.empty())
// //     {
// //         it = graph.erase(it);
// //     }
// //     else
// //     {
// //         it++;
// //     }
// // }
// vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
// // find connected components
// // vector<vector<int>> connectedComponents = FindConnectedComponents(kct);
// vector<pair<Weight, vector<int>>> result;
// // vector<Weight> result_score;

// int ccId = 0;
// for (auto &C : connectedComponents)
// {
//     // cout << "Connected Component: " << ccId++ << endl;
//     // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
//     vector<int> M;
//     // double M_score = numeric_limits<double>::max();
//     /*for compare with another method, there fixed the order of traversed  vertices.*/
//     sort(C.begin(), C.end());
//     NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim);
// }
// printf_results(result, seq2att);
// // auto end = std::chrono::high_resolution_clock::now();
// // std::chrono::duration<double> duration = end - start;
// // cout << "search time: " << duration.count() << "s" << endl;

// // for (int i = 0; i < result.size(); i++)
// // {
// //     cout << "Weight: " << result[i].first.toString() << ", ";
// //     for (int v : result[i].second)
// //     {
// //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
// //         cout << v << ": " << seq2att[v].first << " ";
// //     }
// //     cout << endl;
// //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
// // }
// delete td;
// delete ks;
// return move(result);
// }
bool loadIndexByKey(int k, const std::string &basePath,
                    map<Weight, vector<Edge>, greater<Weight>> &sim_edges,
                    unordered_map<int, vector<int>> &att_vertex,
                    unordered_map<int, vector<int>> &ver_att)
{
    int k_value = 0;
    ifstream k_file(basePath + "k.txt");
    if (!k_file)
    {
        cout << "Cannot open k. " << endl;
        return false;
    }
    vector<int> k_values;
    string line;
    while (getline(k_file, line))
    {
        istringstream iss(line);
        int val;
        while (iss >> val)
        {
            k_values.push_back(val);
        }
    }
    k_file.close();
    for (auto it = k_values.rbegin(); it != k_values.rend(); ++it)
    {
        if (*it >= k)
        {
            k_value = *it;
            break;
        }
    }
    if (k_value == 0)
    {

        return false;
    }

    ifstream sim_file(basePath + to_string(k_value) + "_att.txt");
    if (!sim_file)
    {
        cout << "Cannot open attribute file: " << k_value << endl;
        return false;
    }

    // string line;
    Weight current_weight;
    while (getline(sim_file, line))
    {
        if (line.substr(0, 7) == "Weight:")
        {
            // 假设 Weight 可以由字符串构造
            string weight_str = line.substr(8);
            string numerator = weight_str.substr(0, weight_str.find("/"));
            string denominator = weight_str.substr(weight_str.find("/") + 1, weight_str.find("\n") - weight_str.find("/") - 1);
            current_weight = Weight(stoi(numerator), stoi(denominator)); // 根据实际 Weight 构造函数调整
            sim_edges[current_weight] = vector<Edge>();
        }
        else
        {
            istringstream iss(line);
            int u, v;
            if (iss >> u >> v)
            {
                sim_edges[current_weight].emplace_back(u, v);
            }
        }
    }
    sim_file.close();

    ifstream att_file(basePath + "_att_index_" + to_string(k_value) + ".txt");
    if (!att_file)
    {
        cout << "Cannot open att_index: " << k_value << endl;
        return false;
    }

    int current_subkey = -1;
    while (getline(att_file, line))
    {
        if (line.substr(0, 8) == "SubKey: ")
        {
            current_subkey = stoi(line.substr(8));
            att_vertex[current_subkey] = vector<int>();
        }
        else
        {
            istringstream iss(line);
            int val;
            while (iss >> val)
            {
                att_vertex[current_subkey].push_back(val);
            }
        }
    }
    att_file.close();

    ifstream ver_file(basePath + "_ver_index_" + to_string(k_value) + ".txt");
    if (!ver_file)
    {
        cout << "Cannot open ver_index: " << k_value << endl;
        return false;
    }

    current_subkey = -1;
    while (getline(ver_file, line))
    {
        if (line.substr(0, 8) == "SubKey: ")
        {
            current_subkey = stoi(line.substr(8));
            ver_att[current_subkey] = vector<int>();
        }
        else
        {
            istringstream iss(line);
            int val;
            while (iss >> val)
            {
                ver_att[current_subkey].push_back(val);
            }
        }
    }
    ver_file.close();

    return true;
}
DataGraph LoadGraph(string StructFileName, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &id2att, unordered_map<int, unordered_map<int, Weight>> &s_edges)
{
    DataGraph datagraph;
    ifstream sin(StructFileName.c_str());

    if (!sin)
    {
        cout << "Fail to read " << StructFileName << "." << endl;
        return datagraph;
    }
    string sline;
    while (getline(sin, sline)) // 默认数据集中边不重复
    {
        if (sline.find('#') != string::npos)
            continue;
        std::istringstream lineStream(sline);
        std::vector<std::string> tokens;
        string token;
        while (std::getline(lineStream, token, '\t'))
        {
            tokens.push_back(token);
        }
        int src = stoi(tokens[0]); // string2float
        int dst = stoi(tokens[1]);
        Weight fraction;
        string &fractionStr = tokens[2];
        if (fractionStr == "0/1")
            continue;
        if (id2att.count(src) && id2att.count(dst))
        {
            datagraph.addEdgeNoMatinC(src, dst);

            size_t slashPos = fractionStr.find('/');
            if (slashPos != std::string::npos)
            {
                fraction.numerator = std::stoi(fractionStr.substr(0, slashPos));
                fraction.denominator = std::stoi(fractionStr.substr(slashPos + 1));
            }
            else
            {
                // 如果不是分数格式，直接作为整数处理
                fraction.numerator = std::stoi(fractionStr);
                fraction.denominator = 1;
            }
            int src_id = datagraph.id2seq[src];
            int dst_id = datagraph.id2seq[dst];
            if (src_id > dst_id)
            {
                swap(src_id, dst_id);
            }
            s_edges[src_id][dst_id] = fraction;
        }
    }

    cout << "Loaded dataset successfully!" << endl;
    return move(datagraph);
}
// 索引的，按得分增量扩大，对每个阶段进行枚举
map<Weight, vector<vector<int>>, greater<Weight>> Index_Enumerate(int methods, string dataset, vector<int> posKws, vector<int> negKws,
                                                                  int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    auto greedy_start_time = std::chrono::high_resolution_clock::now();
    Weight opt = WEIGHT_ZERO;
    map<Weight, vector<vector<int>>, greater<Weight>> result;
    map<Weight, vector<Edge>, greater<Weight>> sim_edges;
    unordered_map<int, vector<int>> att_vertex;
    unordered_map<int, vector<int>> ver_att;
    string base_path = "DataGraph/" + dataset + "/index/";
    string StructFileName = base_path + "k.txt";
    // string StructFileName = base_path + to_string(k) + ".txt"; // src, dst, kscore
    // ifstream sin(StructFileName.c_str());
    // // auto start_time = std::chrono::high_resolution_clock::now();

    // if (!sin)
    // {
    //     cout << "Fail to read " << StructFileName << ". Creating..." << endl;
    //     // truss decomposition
    //     // read graph

    //     // construct index
    //     // map<int, shared_ptr<ETNode>> kScore = AdvancedIndex(dg);
    //     Index index(dataset, base_path);
    //     auto end = std::chrono::high_resolution_clock::now();
    //     std::chrono::duration<double> duration = end - start_time;
    //     double search_time = duration.count();
    //     cout << "Construct index finish. Cost " << search_time << "s" << endl;
    //     // index.saveIndexToFile(StructFileName);
    //     // map<Weight, vector<Edge>, greater<Weight>> sim_edges = index.sim_index[k];
    //     // unordered_map<int, vector<int>> att_vertex = index.att_index[k];
    //     // unordered_map<int, vector<int>> ver_att = index.ver_index[k];
    // }
    // sin.close();
    // read k.txt
    int k_value = 0;
    ifstream k_file(StructFileName);
    if (!k_file)
    {
        cout << "Cannot open k. " << endl;
        return result;
    }
    vector<int> k_values;
    string line;

    if (std::getline(k_file, line))
    {
        std::istringstream ss(line);
        std::string token;

        while (std::getline(ss, token, ','))
        {
            try
            {
                k_values.push_back(std::stoi(token));
            }
            catch (const std::invalid_argument &e)
            {
                std::cerr << "Invalid argument: " << e.what() << std::endl;
            }
        }
    }
    k_file.close();
    for (auto it = k_values.rbegin(); it != k_values.rend(); ++it)
    {
        if (*it >= k)
        {
            k_value = *it;
            break;
        }
    }
    if (k_value == 0)
        return result;
    // auto start_time = std::chrono::high_resolution_clock::now();
    string att_text = base_path + to_string(k_value) + "_att.txt";
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att = LoadAtt(0, att_text, posKws, negKws);
    string graph_text = base_path + to_string(k_value) + ".txt";
    unordered_map<int, unordered_map<int, Weight>> s_edges;
    DataGraph dg = LoadGraph(graph_text, id2att, s_edges);
    chrono::duration<double> duration = chrono::high_resolution_clock::now() - start_time;
    cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    // find connected ktruss
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    duration = std::chrono::high_resolution_clock::now() - start_time;
    // cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return result;
    }
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();

    map<Weight, vector<pair<int, int>>, CompareKeys> s_truss_edges;
    for (auto &it : kct)
    {
        int u = it.first;
        for (auto it2 : kct[u])
        {
            int v = it2.first;
            if (u > v)
                continue;

            // if (u_v == WEIGHT_ZERO)
            //     continue;
            Weight w = s_edges[u][v];
            if (!s_truss_edges.count(w))
            {
                s_truss_edges.emplace(w, vector<pair<int, int>>());
                // score_edges.emplace(u_v, unordered_map<int, set<int>>());
            }
            s_truss_edges[w].push_back(pair<int, int>(u, v));
        }
    }
    // map<double, unordered_map<int, set<int>>, CompareKeys> score_edges;

    // cout << "GNKCS start" << endl;
    // Weight bottom_score = GNKCS(methods != 0, true, s_truss_edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    unordered_map<int, unordered_map<int, Weight>> graph_sim;
    Weight bottom_score = GNKCS(true, false, s_truss_edges, seq2att, k, sthd1, sthd2, graph_sim, start_time);
    // duration = std::chrono::high_resolution_clock::now() - start_time;
    // // cout << "Time used: " << duration.count() << "s" << endl;
    // if (duration.count() > max_time)
    // {
    //     return result;
    // }
    // cout << "sim vertex: " << graph_sim.size() << endl;
    // unordered_map<int, unordered_map<int, int>> graph;
    // unordered_map<int, unordered_map<int, Weight>> graph_sim;
    // Weight bottom_score = GNKCS(false, score_edges, seq2att, k, sthd1, sthd2, graph_sim);

    cout << "bottom_score: " << bottom_score.numerator << " / " << bottom_score.denominator << endl;
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration1 = end1 - greedy_start_time;
    cout << "greedy time: " << duration1.count() << "s" << endl;
    duration1 = end1 - start_time;
    if (duration1.count() > max_time)
    {
        return result;
    }
    opt = max(opt, bottom_score);
    unordered_map<int, unordered_map<int, int>> graph;
    for (auto &it_se : s_truss_edges)
    {
        if (it_se.first < opt)
        {
            break;
        }

        for (auto &it1 : it_se.second)
        {
            update_sup(graph, it1);
        }
    }
    vector<Edge> temp_sup_less_edges;
    for (auto &it : graph)
    {
        int u = it.first;
        for (auto &it1 : it.second)
        {
            int v = it1.first;
            if (v > u)
                continue;
            int sup = it1.second;
            if (sup < k - 2)
            {
                temp_sup_less_edges.push_back(pair<int, int>(u, v));
            }
        }
    }
    MaintainKTruss_DeleteEdge(graph, temp_sup_less_edges, k);
    // printGraphAsEdgeList(graph);
    for (auto it = graph.begin(); it != graph.end();)
    {
        // 检查删除边后是否产生孤立点
        if (it->second.empty())
        {
            it = graph.erase(it);
        }
        else
        {
            it++;
        }
    }
    vector<vector<int>> connectedComponents = FindConnectedComponents(graph);
    // find connected components
    // vector<vector<int>> connectedComponents = FindConnectedComponents(kct);

    // vector<Weight> result_score;

    int ccId = 0;
    for (auto &C : connectedComponents)
    {
        for(int u : C){
            cout<< dg.seq2id[u] << ",";
        }
        cout << endl;
        // cout << "Connected Component: " << ccId++ << endl;
        // enumerate all subgraphs and judge whether it meets the nsize and k, calculate the kscore.
        vector<int> M;
        // double M_score = numeric_limits<double>::max();
        /*for compare with another method, there fixed the order of traversed  vertices.*/
        sort(C.begin(), C.end());
        NaiveEnum_AttCntTrussNsizeSortNotConst(WEIGHT_MAX, M, C, graph, seq2att, k, opt, sthd1, sthd2, result, graph_sim, start_time);
    }
    printf_results(dg, result, seq2att);
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;

    // for (int i = 0; i < result.size(); i++)
    // {
    //     cout << "Weight: " << result[i].first.toString() << ", ";
    //     for (int v : result[i].second)
    //     {
    //         // cout << dg.seq2id[v] << ": " << seq2att[v].first << " ";
    //         cout << v << ": " << seq2att[v].first << " ";
    //     }
    //     cout << endl;
    //     // cout << "score: " << result_score[i].numerator << '\\' << result_score[i].denominator << endl;
    // }
    delete td;
    return move(result);
}
