#include "KC.h"
#include <climits>
bool IsAllAttrIn(unordered_map<int, unordered_map<int, int>> &G, vector<int> &Kws, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    for (auto attr : Kws)
    {
        bool isAttrIn = false;
        for (auto node : G)
        {
            int node_id = node.first;
            vector<int> att = seq2att.at(node_id).second.second;
            if (count(att.begin(), att.end(), attr))
            {
                isAttrIn = true;
                break;
            }
        }
        if (!isAttrIn)
        {
            return false;
        }
    }
    return true;
}
int GetKDistU(unordered_map<int, unordered_map<int, int>> &G, int u, vector<int> &queryAttrs, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    unordered_set<int> queryAttrsH(queryAttrs.begin(), queryAttrs.end());
    queue<pair<int, int>> Q; // <node, distance>
    unordered_set<int> visited;

    Q.push(make_pair(u, 0));
    visited.emplace(u);

    while (!Q.empty())
    {
        int node = Q.front().first;
        int dist = Q.front().second;
        Q.pop();

        // 检查节点属性
        auto &nodeAttrs = seq2att.at(node).second.second;
        for (auto attr : nodeAttrs)
        {
            queryAttrsH.erase(attr);
        }

        // 如果所有查询属性都已找到
        if (queryAttrsH.empty())
        {
            return dist;
        }

        // 遍历邻居节点
        for (auto &nbr_entry : G.at(node))
        {
            int nbr = nbr_entry.first;
            if (visited.find(nbr) == visited.end())
            {
                Q.push(make_pair(nbr, dist + 1));
                visited.emplace(nbr);
            }
        }
    }

    return INT_MAX;
}
unordered_map<int, unordered_map<int, int>> GetMaxKTruss(
    unordered_map<int, unordered_map<int, int>> &graph,
    int k)
{
    vector<Edge> supUqEdges;

    // 计算新添加边的支持度
    int supQEdgeCnt = 0;
    for (const auto &u_v : graph)
    {
        int u = u_v.first;
        for (auto &v_sup : graph[u])
        {

            int v = v_sup.first;
            if (u < v)
                continue;
            int support = 0;

            // 确保u的度数小于等于v的度数
            // if (graph[u].size() > graph[v].size())
            //     swap(u, v);

            // 计算共同邻居数量
            for (const auto &nbr_entry : graph[u])
            {
                int nbr = nbr_entry.first;
                if (graph[v].count(nbr))
                {
                    support++;
                }
            }

            graph[u][v] = support;
            graph[v][u] = support;

            if (support < k - 2)
            {

                Edge e = GetEdge(u, v);
                supUqEdges.push_back(e);
            }
            else
            {
                supQEdgeCnt++;
            }
        }
    }

    // 处理边集
    size_t i = 0;
    for (; i < supUqEdges.size(); i++)
    {
        if (supQEdgeCnt == 0)
        {
            return {};
        }

        Edge e = supUqEdges[i];
        int u = e.first, v = e.second;

        // 确保u的度数小于等于v的度数
        if (graph[u].size() > graph[v].size())
            swap(u, v);

        // 遍历u的邻居
        for (const auto &nbr_entry : graph[u])
        {
            int nbr = nbr_entry.first;
            if (!graph[v].count(nbr))
            {
                continue;
            }

            Edge nbr_v = GetEdge(nbr, v);
            Edge nbr_u = GetEdge(nbr, u);

            // 更新nbr_v边的支持度
            if (graph[nbr][v] != -1)
            {
                int &support = graph[nbr][v];
                if (support == k - 2)
                {
                    supQEdgeCnt--;
                    supUqEdges.push_back(nbr_v);
                }
                support--;
                graph[v][nbr] = support;
            }

            // 更新nbr_u边的支持度
            if (graph[nbr][u] != -1)
            {
                int &support = graph[nbr][u];
                if (support == k - 2)
                {
                    supQEdgeCnt--;
                    supUqEdges.push_back(nbr_u);
                }
                support--;
                graph[u][nbr] = support;
            }
        }

        // 从图中删除这条边
        graph[u].erase(v);
        graph[v].erase(u);
    }

    // 检查weightEqPairs是否为空
    // bool weightPairsEmpty = true;
    // for (const auto &u_entry : weightEqPairs)
    // {
    //     if (!u_entry.second.empty())
    //     {
    //         weightPairsEmpty = false;
    //         break;
    //     }
    // }
    // if (weightPairsEmpty)
    // {
    //     return {};
    // }

    // 删除孤立节点
    vector<int> delNodes;
    for (const auto &u_entry : graph)
    {
        if (u_entry.second.empty())
        {
            delNodes.push_back(u_entry.first);
        }
    }

    for (int node : delNodes)
    {
        graph.erase(node);
    }

    return move(graph);
}
map<Weight, vector<vector<int>>, greater<Weight>> KC(string dataset, vector<int> posKws,
                                                     int k, chrono::high_resolution_clock::time_point start_time)
{
    // auto start_time = std::chrono::high_resolution_clock::now();
    map<Weight, vector<vector<int>>, greater<Weight>> KCresults;
    unordered_map<int, unordered_map<int, int>> KCResult;
    // auto start = std::chrono::high_resolution_clock::now();
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att;
    int resultCnt = 0;
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    id2att = LoadAtt(1, AttFileName, posKws, vector<int>());

    DataGraph dg = LoadGraph(dataset, id2att);
    // find connected ktruss
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    if (kct.empty())
    {
        cout << "No k-truss found!" << endl;
        return KCresults;
    }
    else
    {
        cout << "k-truss found!" << endl;
    }
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    map<Weight, vector<pair<int, int>>, CompareKeys> score_edges;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();
    if (!IsAllAttrIn(kct, posKws, seq2att))
    {
        return KCresults;
    }
    // 计算每个节点的k距离
    int d = -1;
    unordered_map<int, int> kDistMap;
    for (auto &node_entry : kct)
    {
        int node = node_entry.first;
        int dist = GetKDistU(kct, node, posKws, seq2att);
        kDistMap.emplace(node, dist);
        if (dist > d)
        {
            d = dist;
        }
    }
    KCResult = kct;
    while (IsAllAttrIn(kct, posKws, seq2att))
    {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - start_time;
        if (duration.count() > max_time)
        {
            return {};
        }
        // 删除k距离大于等于d的节点
        vector<int> delNodes;
        for (auto &node_entry : kct)
        {
            int node = node_entry.first;
            if (kDistMap.at(node) >= d)
            {
                delNodes.push_back(node);
            }
        }

        for (auto node : delNodes)
        {
            for (auto &neighbor_entry : kct[node])
            {
                kct[neighbor_entry.first].erase(node);
            }
            kct.erase(node);
        }

        // 重新计算最大k-truss
        kct = GetMaxKTruss(kct, k);
        if (kct.empty())
        {
            break;
        }

        // 更新k距离
        int d1 = -1;
        for (auto &node_entry : kct)
        {
            int node = node_entry.first;
            int dist = GetKDistU(kct, node, posKws, seq2att);
            kDistMap[node] = dist;
            if (dist > d1)
            {
                d1 = dist;
            }
        }

        // 更新结果
        if (d1 < d)
        {
            d = d1;
            KCResult = kct;
        }
    }

    // 构建结果
    vector<vector<int>> results = FindConnectedComponents(KCResult);
    if(results.empty()){
        return {};
    }
    for (auto &result : results)
    {
        vector<int> idResult;
        idResult.reserve(result.size());

        for (int seq : result)
        {
            int id = dg.seq2id[seq];
            idResult.push_back(id);
        }

        KCresults[Weight(d, 1)].push_back(move(idResult));
    }
    return move(KCresults);
}