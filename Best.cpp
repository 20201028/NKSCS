#include "Best.h"
#include <climits>
void MaintainSupport_AddEdge(unordered_map<int, unordered_map<int, int>> &graph, vector<set<Edge>> &nonDecSup, Edge e, int k)
{
    int u = e.first, v = e.second;

    // 检查节点是否存在
    bool u_exists = graph.count(u) > 0;
    bool v_exists = graph.count(v) > 0;

    if (u_exists && v_exists)
    {
        // 确保u的度数小于等于v的度数
        if (graph[u].size() > graph[v].size())
            swap(u, v);

        // 遍历u的所有邻居
        int com = 0;
        for (auto &nbr_entry : graph[u])
        {
            int nbr = nbr_entry.first;
            // 检查nbr是否是v的邻居
            if (graph[v].count(nbr) > 0)
            {
                com++;
                Edge nbr_v = GetEdge(nbr, v);
                if (++graph[nbr][v] == k - 2)
                {
                    nonDecSup[0].erase(nbr_v);
                    nonDecSup[1].insert(nbr_v);
                }
                ++graph[v][nbr];
                Edge nbr_u = GetEdge(nbr, u);
                if (++graph[nbr][u] == k - 2)
                {
                    nonDecSup[0].erase(nbr_u);
                    nonDecSup[1].insert(nbr_u);
                }
                ++graph[u][nbr];
            }
        }

        graph[u][v] = com;
        graph[v][u] = com;
        if (com >= k - 2)
        {
            nonDecSup[1].insert(e);
        }
        else
        {
            nonDecSup[0].insert(e);
        }
    }
    else
    {
        // 添加边到图中，支持度初始为0
        graph[u][v] = 0;
        graph[v][u] = 0;
        nonDecSup[0].insert(e);
    }
}
set<Edge> GetMaxKTrussInc(
    unordered_map<int, unordered_map<int, int>> &graph,
    int k,
    vector<set<Edge>> &nonDecSup)
{
    set<Edge> result;
    if (nonDecSup[1].empty())
    {
        return result;
    }

    // 收集所有符合条件的边
    vector<Edge> newEdges, supUqEdges;
    newEdges.reserve(nonDecSup[1].size());
    supUqEdges.reserve(nonDecSup[1].size());

    // 将这些边添加到图中
    for (const auto &e : nonDecSup[1])
    {
        newEdges.push_back(e);
        graph[e.first][e.second] = -1; // 初始化支持度为0
        graph[e.second][e.first] = -1;
    }

    // 计算新添加边的支持度
    int supQEdgeCnt = 0;
    for (const auto &e : newEdges)
    {
        int u = e.first, v = e.second;
        int support = 0;

        // 确保u的度数小于等于v的度数
        if (graph[u].size() > graph[v].size())
            swap(u, v);

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

            supUqEdges.push_back(e);
        }
        else
        {
            supQEdgeCnt++;
        }
    }

    // 处理边集
    for (size_t i = 0; i < supUqEdges.size(); i++)
    {
        if (supQEdgeCnt == 0)
        {
            break;
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
        // if(graph[u].size() == 0)graph.erase(u);
        graph[v].erase(u);
        // if(graph[v].size() == 0)graph.erase(v);
    }

    // 收集结果

    for (const auto &e : newEdges)
    {
        if (graph[e.first].count(e.second))
        {
            result.insert(e);
            graph[e.second][e.first] = -1;
            graph[e.first][e.second] = -1;
            nonDecSup[1].erase(e);
        }
    }

    return move(result);
}
void AddEdges(set<Edge> &edges, unordered_map<int, unordered_map<int, int>> &T,
              vector<unordered_map<int, unordered_map<int, int>>> &maxCTrusses,
              map<Edge, int> &CTrussOfEdge, map<int, vector<int>> &newPairs,
              queue<int> &emptyIndex)
{
    int index = maxCTrusses.size();
    unordered_map<int, unordered_map<int, int>> TCTrussGraph;
    vector<vector<Edge>> newTCs;

    // do BFS, when traversing to existing edges, stop further traverse
    unordered_set<int> mergeTCTrusses;
    set<Edge> visited;

    while (!edges.empty())
    {
        vector<Edge> result;
        queue<Edge> Q;
        Edge beginE = *edges.begin();
        Q.push(beginE);
        visited.insert(beginE);

        while (!Q.empty())
        {
            Edge e = Q.front();
            Q.pop();

            auto it = CTrussOfEdge.find(e);
            if (it != CTrussOfEdge.end())
            {
                mergeTCTrusses.insert(it->second);
                continue;
            }

            int u = e.first, v = e.second;
            if (T[u].size() > T[v].size())
            {
                swap(u, v);
            }

            result.push_back(e);
            edges.erase(e);

            // triangle connected Find common neighbors
            // vector<int> comNeighbors;
            // for (const auto &nbr_entry : T[u])
            // {
            //     int nbr = nbr_entry.first;
            //     if (T[v].count(nbr))
            //     {
            //         comNeighbors.push_back(nbr);
            //     }
            // }
            // for (auto w : comNeighbors)
            // {
            //     Edge w_v = GetEdge(w, v);
            //     if (visited.insert(w_v).second)
            //     {
            //         Q.push(w_v);
            //     }
            //     Edge w_u = GetEdge(w, u);
            //     if (visited.insert(w_u).second)
            //     {
            //         Q.push(w_u);
            //     }
            // }
            for (const auto &nbr_entry : T[u])
            {
                int w = nbr_entry.first;
                Edge w_v = GetEdge(u, w);
                if (visited.insert(w_v).second)
                {
                    Q.push(w_v);
                }
            }
            for (const auto &nbr_entry : T[v])
            {
                int w = nbr_entry.first;
                Edge w_v = GetEdge(v, w);
                if (visited.insert(w_v).second)
                {
                    Q.push(w_v);
                }
            }
        }

        newTCs.push_back(result);
        TCTrussGraph[index][index] = 0; // Add self-edge

        for (int TCIndex : mergeTCTrusses)
        {
            TCTrussGraph[index][TCIndex] = 1;
            TCTrussGraph[TCIndex][index] = 1;
        }
        index++;

        mergeTCTrusses.clear();
        visited.clear();
    }

    // Find connected components in TCTrussGraph
    vector<vector<int>> newTCIndexes;
    unordered_set<int> visitedNodes;

    for (const auto &node_entry : TCTrussGraph)
    {
        int node = node_entry.first;
        if (visitedNodes.count(node))
            continue;

        vector<int> component;
        queue<int> q;
        q.push(node);
        visitedNodes.insert(node);

        while (!q.empty())
        {
            int current = q.front();
            q.pop();
            component.push_back(current);

            for (const auto &neighbor_entry : TCTrussGraph[current])
            {
                int neighbor = neighbor_entry.first;
                if (!visitedNodes.count(neighbor))
                {
                    visitedNodes.insert(neighbor);
                    q.push(neighbor);
                }
            }
        }
        sort(component.begin(), component.end());
        newTCIndexes.push_back(component);
    }

    // Algorithm 5
    int TCTrussSize = maxCTrusses.size();
    for (auto &newTC : newTCIndexes)
    {
        unordered_map<int, unordered_map<int, int>> result; // T_N
        vector<int> newNodes;                               // N
        vector<vector<int>> preTCTrussNodes;

        for (auto index : newTC)
        {
            if (index < TCTrussSize)
            {
                auto &graph = maxCTrusses[index];
                for (const auto &u_entry : graph)
                {
                    int u = u_entry.first;
                    for (const auto &v_entry : u_entry.second)
                    {
                        int v = v_entry.first;
                        result[u][v] = -1; // Default edge value
                    }
                }

                vector<int> nodes;
                for (const auto &u_entry : graph)
                {
                    nodes.push_back(u_entry.first);
                }
                preTCTrussNodes.push_back(nodes);
                emptyIndex.push(index);
            }
            else
            {
                for (auto &edge : newTCs[index - TCTrussSize])
                {
                    int u = edge.first, v = edge.second;
                    if (result.find(u) == result.end())
                    {
                        newNodes.push_back(u);
                    }
                    if (result.find(v) == result.end())
                    {
                        newNodes.push_back(v);
                    }
                    result[u][v] = -1;
                    result[v][u] = -1;
                }
            }
        }

        // AddTCTruss(result);
        int index = emptyIndex.empty() ? maxCTrusses.size() : emptyIndex.front();

        if (!emptyIndex.empty())
        {
            emptyIndex.pop();
            maxCTrusses[index] = result;
        }
        else
        {
            maxCTrusses.push_back(result);
        }

        // 遍历新TC-truss中的所有边
        for (const auto &u_entry : result)
        {
            int u = u_entry.first;
            for (const auto &v_entry : u_entry.second)
            {
                int v = v_entry.first;
                // 确保每条边只处理一次（避免重复处理u-v和v-u）
                if (u < v)
                {
                    CTrussOfEdge[make_pair(u, v)] = index;
                }
            }
        }

        std::unordered_set<int> oldNodeSet_1;
        std::unordered_map<int, int> nodeCntMap;

        for (const auto &nodes : preTCTrussNodes)
        {
            for (int node : nodes)
            {
                oldNodeSet_1.insert(node);
                nodeCntMap[node]++;
            }
        }

        int newNodeCnt = newNodes.size();
        for (int i = 0; i < newNodeCnt; ++i)
        {
            for (int j = i + 1; j < newNodeCnt; ++j)
            {
                int u = newNodes[i], v = newNodes[j];
                if (u < v)
                {
                    newPairs[u].push_back(v);
                }
                else
                {
                    newPairs[v].push_back(u);
                }
            }
        }

        for (int i = 0; i < newNodeCnt; ++i)
        {
            int u = newNodes[i];
            for (int v : oldNodeSet_1)
            {
                if (u < v)
                {
                    newPairs[u].push_back(v);
                }
                else
                {
                    newPairs[v].push_back(u);
                }
            }
        }

        for (const auto &oldNodes : preTCTrussNodes)
        {
            std::vector<int> diffNodes;
            std::unordered_set<int> oldNodeSet_2(oldNodeSet_1.begin(), oldNodeSet_1.end());

            for (int oldNode : oldNodes)
            {
                nodeCntMap[oldNode]--;
                if (nodeCntMap[oldNode] == 0)
                {
                    diffNodes.push_back(oldNode);
                    oldNodeSet_1.erase(oldNode);
                }
                oldNodeSet_2.erase(oldNode);
            }

            for (int u : diffNodes)
            {
                for (int v : oldNodeSet_2)
                {
                    if (u < v)
                    {
                        newPairs[u].push_back(v);
                    }
                    else
                    {
                        newPairs[v].push_back(u);
                    }
                }
            }
        }
    }
}
// void AddEdges(set<Edge> &edges, unordered_map<int, unordered_map<int, int>> &T,
//               vector<unordered_map<int, unordered_map<int, int>>> &maxCTrusses,
//               map<Edge, int> &CTrussOfEdge, map<int, vector<int>> &newPairs,
//               queue<int> &emptyIndex)
// {
//     int index = maxCTrusses.size();
//     unordered_map<int, unordered_map<int, int>> CTrussGraph;
//     vector<vector<Edge>> newCs;

//     // do BFS, when traversing to existing edges, stop further traverse
//     unordered_set<int> mergeCTrusses;
//     set<Edge> visited;
//     while (!edges.empty())
//     {
//         vector<Edge> result;
//         queue<Edge> Q;
//         Edge beginE = *edges.begin();
//         Q.push(beginE);
//         visited.insert(beginE);
//         while (!Q.empty())
//         {
//             Edge e = Q.front();
//             Q.pop();

//             auto it = CTrussOfEdge.find(e);
//             if (it != CTrussOfEdge.end())
//             {
//                 mergeCTrusses.insert(it->second);
//                 continue;
//             }

//             int u = e.first, v = e.second;
//             if (T[u].size() > T[v].size())
//             {
//                 swap(u, v);
//             }

//             result.push_back(e);
//             edges.erase(e);

//             for (auto w : T[u])
//             {
//                 Edge w_u = GetEdge(w.first, u);
//                 if (visited.insert(w_u).second)
//                 {
//                     Q.push(w_u);
//                 }
//             }
//             for (auto w : T[v])
//             {
//                 Edge w_v = GetEdge(w.first, v);
//                 if (visited.insert(w_v).second)
//                 {
//                     Q.push(w_v);
//                 }
//             }
//         }

//         newCs.push_back(result);
//         CTrussGraph.emplace(index, unordered_map<int, int>());
//         for (int TCIndex : mergeCTrusses)
//         {
//             if (!CTrussGraph.count(TCIndex))
//             {
//                 CTrussGraph.emplace(TCIndex, unordered_map<int, int>());
//             }
//             CTrussGraph[index].emplace(TCIndex, -1);
//             CTrussGraph[TCIndex].emplace(index, -1);
//         }
//         index++;

//         mergeCTrusses.clear();
//         visited.clear();
//     }

//     vector<vector<int>> newCIndexes;
//     // set<int> VisitedNId;
//     unordered_set<int> visitedNodes;
//     for (auto it : CTrussGraph)
//     {
//         if (it.second.size() == 0)
//         {
//             int NId = it.first;
//             visitedNodes.insert(NId);
//             vector<int> component;
//             component.push_back(NId);
//             newCIndexes.push_back(component);
//         }
//     }
//     for (auto &node_entry : CTrussGraph)
//     {
//         int node = node_entry.first;
//         if (visitedNodes.count(node))
//             continue;

//         vector<int> component;
//         queue<int> q;
//         q.push(node);
//         visitedNodes.insert(node);

//         while (!q.empty())
//         {
//             int current = q.front();
//             q.pop();
//             component.push_back(current);

//             for (auto &neighbor_entry : CTrussGraph[current])
//             {
//                 int neighbor = neighbor_entry.first;
//                 if (!visitedNodes.count(neighbor))
//                 {
//                     visitedNodes.insert(neighbor);
//                     q.push(neighbor);
//                 }
//             }
//         }

//         sort(component.begin(), component.end());
//         newCIndexes.push_back(component);
//     }
//     // vector<int> CcNIdV(1, 0);
//     // for (auto it : CTrussGraph)
//     // {
//     //     if (it.second.size() == 0)
//     //     {
//     //         int NId = it.first;
//     //         VisitedNId.insert(NId);
//     //         CcNIdV[0] = NId;
//     //         newCIndexes.push_back(CcNIdV);
//     //     }
//     // }

//     // for (auto it : CTrussGraph)
//     // {
//     //     int NId = it.first;
//     //     if (!VisitedNId.count(NId))
//     //     {
//     //         VisitedNId.insert(NId);
//     //         queue<int> NIdQ;
//     //         NIdQ.push(NId);
//     //         CcNIdV.clear();
//     //         CcNIdV.push_back(NId);
//     //         while (!NIdQ.empty())
//     //         {
//     //             int Node = NIdQ.front();
//     //             NIdQ.pop();
//     //             for(auto it : CTrussGraph[Node]){
//     //                 int OutNId = it.first;
//     //                 if (!VisitedNId.count(OutNId))
//     //                 {
//     //                     NIdQ.push(OutNId);
//     //                     VisitedNId.insert(OutNId);
//     //                     CcNIdV.push_back(OutNId);
//     //                 }
//     //             }

//     //         }
//     //         sort(CcNIdV.begin(), CcNIdV.end());
//     //         newCIndexes.push_back(CcNIdV);
//     //     }
//     // }

//     // Algorithm 5

//     int CTrussSize = maxCTrusses.size();
//     for (auto &newC : newCIndexes)
//     {
//         unordered_map<int, unordered_map<int, int>> result; // T_N
//         vector<int> newNodes;                               // N
//         vector<vector<int>> preCTrussNodes;

//         for (auto &idx : newC)
//         {
//             if (idx < CTrussSize)
//             {
//                 // Merge existing C-truss
//                 auto &graph = maxCTrusses[idx];
//                 for (auto &u_entry : graph)
//                 {
//                     int u = u_entry.first;
//                     for (auto &v_entry : u_entry.second)
//                     {
//                         int v = v_entry.first;
//                         result[u][v] = v_entry.second;
//                         result[v][u] = v_entry.second;
//                     }
//                 }

//                 // Collect nodes from this C-truss
//                 vector<int> nodes;
//                 for (auto &u_entry : graph)
//                 {
//                     nodes.push_back(u_entry.first);
//                 }
//                 preCTrussNodes.push_back(nodes);
//             }
//             else
//             {
//                 // Add new edges
//                 for (auto &edge : newCs[idx - CTrussSize])
//                 {
//                     int u = edge.first, v = edge.second;
//                     result[u][v] = 0; // Initialize support to 0
//                     result[v][u] = 0;

//                     if (find(newNodes.begin(), newNodes.end(), u) == newNodes.end())
//                     {
//                         newNodes.push_back(u);
//                     }
//                     if (find(newNodes.begin(), newNodes.end(), v) == newNodes.end())
//                     {
//                         newNodes.push_back(v);
//                     }
//                 }
//             }
//         }
//         // for (auto &newTC : newCIndexes)
//         // {
//         //     unordered_map<int, unordered_map<int, int>> result; // T_N
//         //     vector<int> newNodes;              // N
//         //     vector<vector<int>> preTCTrussNodes;
//         //     int newCsize = newTC.size();
//         //     preTCTrussNodes.reserve(newCsize);
//         //     // we do the mapping of new and old triangle-connected k-trusses when dynamically computing the new triangle-connected k-trusses
//         //     // so we directly obtain the corresponding \mathcal{T}_O
//         //     for (auto &index : newTC)
//         //     {
//         //         if (index < TCTrussSize)
//         //         {
//         //             unordered_map<int, unordered_map<int, int>> graph = maxCTrusses[index];
//         //             for (TUNGraph::TEdgeI EI = graph->BegEI(); EI != graph->EndEI(); EI++)
//         //             {
//         //                 int u = EI.GetSrcNId(), v = EI.GetDstNId();
//         //                 result->AddNodeUnchecked(u);
//         //                 result->AddNodeUnchecked(v);
//         //                 result->AddEdgeUnchecked(u, v);
//         //             }
//         //             emptyIndex.push(index);
//         //             vector<int> nodes;
//         //             nodes.reserve(graph->GetNodes());
//         //             for (TUNGraph::TNodeI NI = graph->BegNI(); NI != graph->EndNI(); NI++)
//         //             {
//         //                 nodes.push_back(NI.GetId());
//         //             }
//         //             preTCTrussNodes.push_back(nodes);
//         //         }
//         //         else
//         //         {
//         //             for (auto &edge : newCs[index - TCTrussSize])
//         //             {
//         //                 int u = edge.first, v = edge.second;
//         //                 if (!result->IsNode(u))
//         //                 {
//         //                     result->AddNode(u);
//         //                     newNodes.push_back(u);
//         //                 }
//         //                 if (!result->IsNode(v))
//         //                 {
//         //                     result->AddNode(v);
//         //                     newNodes.push_back(v);
//         //                 }
//         //                 result->AddEdge(u, v);
//         //             }
//         //         }
//         //     }
//         int index = emptyIndex.empty() ? maxCTrusses.size() : emptyIndex.front();

//         if (!emptyIndex.empty())
//         {
//             emptyIndex.pop();
//             maxCTrusses[index] = result;
//         }
//         else
//         {
//             maxCTrusses.push_back(result);
//         }

//         // 遍历新TC-truss中的所有边
//         for (const auto &u_entry : result)
//         {
//             int u = u_entry.first;
//             for (const auto &v_entry : u_entry.second)
//             {
//                 int v = v_entry.first;
//                 // 确保每条边只处理一次（避免重复处理u-v和v-u）
//                 if (u < v)
//                 {
//                     CTrussOfEdge[make_pair(u, v)] = index;
//                 }
//             }
//         }

//         unordered_set<int> oldNodeSet_1; // O_1
//         unordered_map<int, int> nodeCntMap;
//         for (auto &nodes : preCTrussNodes)
//         {
//             for (auto node : nodes)
//             {
//                 oldNodeSet_1.insert(node);
//                 nodeCntMap[node]++;
//             }
//         }

//         // obtain P(N)
//         int newNodeCnt = newNodes.size();
//         for (int i = 0; i < newNodeCnt; i++)
//         {
//             for (int j = i + 1; j < newNodeCnt; j++)
//             {
//                 int u = newNodes[i], v = newNodes[j];
//                 if (u < v)
//                 {
//                     newPairs[u].push_back(v);
//                 }
//                 else
//                 {
//                     newPairs[v].push_back(u);
//                 }
//             }
//         }

//         // obtain {(u, v)|u �� N, v �� O_1}
//         for (int i = 0; i < newNodeCnt; i++)
//         {
//             int u = newNodes[i];
//             for (auto &v : oldNodeSet_1)
//             {
//                 if (u < v)
//                 {
//                     newPairs[u].push_back(v);
//                 }
//                 else
//                 {
//                     newPairs[v].push_back(u);
//                 }
//             }
//         }

//         // obtain P(O_1) \ \mathcal{P}(\mathcal{T}_O)
//         for (auto &oldNodes : preCTrussNodes)
//         {
//             vector<int> diffNodes; // D
//             diffNodes.reserve(oldNodes.size());
//             unordered_set<int> oldNodeSet_2(oldNodeSet_1); // O_2
//             for (auto &oldNode : oldNodes)
//             {
//                 nodeCntMap[oldNode]--;
//                 if (nodeCntMap[oldNode] == 0)
//                 {
//                     diffNodes.push_back(oldNode);
//                     oldNodeSet_1.erase(oldNode);
//                 }
//                 oldNodeSet_2.erase(oldNode);
//             }

//             for (auto &u : diffNodes)
//             {
//                 for (auto &v : oldNodeSet_2)
//                 {
//                     if (u < v)
//                     {
//                         newPairs[u].push_back(v);
//                     }
//                     else
//                     {
//                         newPairs[v].push_back(u);
//                     }
//                 }
//             }
//         }
//     }
// }
void GetNonAscDeg(const unordered_map<int, unordered_map<int, int>> &G, vector<pair<int, int>> &deg)
{
    deg.reserve(G.size());

    // 遍历图中的所有节点
    for (const auto &node_entry : G)
    {
        int node = node_entry.first;
        int degree = node_entry.second.size(); // 邻居数量即为度数

        // 存储(度数, 节点ID)对
        deg.emplace_back(degree, node);
    }

    // 按度数非升序排序
    sort(deg.begin(), deg.end(), [](const pair<int, int> &a, const pair<int, int> &b)
         {
             return a.first < b.first; // 降序排序
         });
}
void AddClique(vector<vector<int>> &maxCliques,
               map<int, set<int>> &cliquesOfNode,
               set<int> &notEmptyIndex, set<int> &emptyIndex,
               set<int> &resultCliques,
               vector<int> &clique)
{
    int index;

    // 获取可用的索引位置
    if (emptyIndex.empty())
    {
        index = maxCliques.size();
        maxCliques.push_back(clique);
    }
    else
    {
        index = *emptyIndex.begin();
        maxCliques[index] = clique;
        emptyIndex.erase(emptyIndex.begin());
    }

    // 为clique中的每个节点记录所属的clique索引
    for (auto node : clique)
    {
        cliquesOfNode[node].insert(index);
    }

    // 更新索引集合
    notEmptyIndex.insert(index);
    resultCliques.insert(index);
}
void InsertNode(vector<int> &clique, int node)
{
    int size = clique.size();
    if (size == 0)
    {
        clique.push_back(node);
        return;
    }
    clique.push_back(node);
    sort(clique.begin(), clique.end());
}
pair<int, int> GetMinDegNode(const unordered_map<int, vector<int>> &nodeNeighbors,
                             const vector<int> &cand)
{
    int minDeg = INT_MAX;
    int minDegNode = -1;

    for (auto node : cand)
    {
        int degree = nodeNeighbors.at(node).size();
        if (degree < minDeg)
        {
            minDeg = degree;
            minDegNode = node;
        }
    }

    return make_pair(minDegNode, minDeg);
}
bool MaxEval(vector<int> &cand, vector<int> &E,
             const unordered_map<int, unordered_map<int, int>> &G)
{
    int candSize = cand.size();

    // 存储节点及其邻居信息
    unordered_map<int, vector<int>> nodeNeighbors;
    for (auto node : cand)
    {
        if (G.find(node) != G.end())
        {
            vector<int> neighbors;
            for (const auto &nbr_entry : G.at(node))
            {
                neighbors.push_back(nbr_entry.first);
            }
            nodeNeighbors.emplace(node, move(neighbors));
        }
    }

    // 获取度数最小的节点
    auto node_deg = GetMinDegNode(nodeNeighbors, cand);
    int minDeg = node_deg.second;
    int minDegNode = node_deg.first;

    if (minDeg == candSize - 1)
    {
        return true;
    }

    // 获取最小度数节点的邻居
    vector<int> nbrs(nodeNeighbors[minDegNode].begin(),
                     nodeNeighbors[minDegNode].end());

    // 计算差集
    vector<int> diff;
    sort(nbrs.begin(), nbrs.end());
    sort(E.begin(), E.end());
    set_difference(nbrs.begin(), nbrs.end(),
                   E.begin(), E.end(),
                   back_inserter(diff));

    // 检查每个差异节点
    for (auto node : diff)
    {
        if (G.find(node) != G.end())
        {
            int nodeDeg = G.at(node).size();
            if (nodeDeg >= candSize)
            {
                // 获取节点的所有邻居
                vector<int> uNbrs;
                for (const auto &nbr_entry : G.at(node))
                {
                    uNbrs.push_back(nbr_entry.first);
                }

                // 检查候选集是否包含在邻居中
                sort(uNbrs.begin(), uNbrs.end());
                sort(cand.begin(), cand.end());
                if (includes(uNbrs.begin(), uNbrs.end(),
                             cand.begin(), cand.end()))
                {
                    return false;
                }
            }
        }
    }
    return true;
}
bool MaxEval(std::vector<int> &cand, std::vector<int> &E, std::unordered_map<int, std::vector<int>> &nodeInfo)
{
    int candSize = cand.size();
    auto node_deg = GetMinDegNode(nodeInfo, cand);
    int minDeg = node_deg.second, minDegNode = node_deg.first;
    if (minDeg == candSize - 1)
    {
        return true;
    }

    vector<int> diff;
    auto &nbrs = nodeInfo.at(minDegNode);
    diff.reserve(minDeg);
    set_difference(nbrs.begin(), nbrs.end(), E.begin(), E.end(), back_inserter(diff));

    for (auto u : diff)
    {
        auto &uNbrs = nodeInfo.at(u);
        if (uNbrs.size() >= candSize)
        {
            sort(uNbrs.begin(), uNbrs.end());
            sort(cand.begin(), cand.end());
            if (includes(uNbrs.begin(), uNbrs.end(),
                         cand.begin(), cand.end()))
            {
                return false;
            }
        }
    }
    return true;
}
bool MaxEval(vector<int> &cand, set<int> &E, unordered_map<int, vector<int>> &nodeInfo)
{
    int candSize = cand.size();
    auto node_deg = GetMinDegNode(nodeInfo, cand);
    int minDeg = node_deg.second, minDegNode = node_deg.first;
    if (minDeg == candSize - 1)
    {
        return true;
    }

    auto &nbrs = nodeInfo.at(minDegNode);
    for (auto u : nbrs)
    {
        if (E.find(u) == E.end()) // u �� N(minDegNode) / E
        {
            auto &uNbrs = nodeInfo.at(u);
            if (uNbrs.size() >= candSize)
            {
                sort(uNbrs.begin(), uNbrs.end());
                sort(cand.begin(), cand.end());
                if (
                    includes(uNbrs.begin(), uNbrs.end(),
                             cand.begin(), cand.end()))
                {
                    return false;
                }
            }
        }
    }
    return true;
}
void DelClique(vector<vector<int>> &maxCliques,
               map<int, set<int>> &cliquesOfNode,
               set<int> &notEmptyIndex, set<int> &emptyIndex,
               set<int> &resultCliques, int index)
{
    for (auto node : maxCliques[index])
    {
        cliquesOfNode[node].erase(index);
    }
    notEmptyIndex.erase(index);
    resultCliques.erase(index);
    emptyIndex.insert(index);
    maxCliques[index].clear();
}
void AddEdge(unordered_map<int, unordered_map<int, int>> &G,
             vector<vector<int>> &maxCliques,
             map<int, set<int>> &cliquesOfNode,
             set<int> &notEmptyIndex, set<int> &emptyIndex,
             set<int> &resultCliques,
             int u, int v)
{
    // 直接添加新边 {u, v}
    if (G.find(u) == G.end() || G.find(v) == G.end())
    {
        // 添加边到图中（默认边信息为-1）
        G[u][v] = -1;
        G[v][u] = -1;

        // 创建新的clique {u, v}
        vector<int> clique;
        if (u < v)
        {
            clique = {u, v};
        }
        else
        {
            clique = {v, u};
        }
        AddClique(maxCliques, cliquesOfNode,
                  notEmptyIndex, emptyIndex, resultCliques, clique);
        return;
    }

    // 添加边到图中（默认边信息为-1）
    G[u][v] = -1;
    G[v][u] = -1;

    // 确保d(u) <= d(v)
    if (G[u].size() > G[v].size())
    {
        swap(u, v);
    }

    // 获取v的邻居集合
    set<int> vNbrs;
    for (const auto &nbr_entry : G[v])
    {
        vNbrs.insert(nbr_entry.first);
    }

    // 处理u所在的cliques
    set<int> cliquesU = cliquesOfNode[u];
    set<int> cliquesV = cliquesOfNode[v];
    set<vector<int>> candCliques;

    for (auto cliqueIndex : cliquesU)
    {
        auto &clique = maxCliques[cliqueIndex];
        vector<int> cand;
        cand.reserve(min(clique.size(), vNbrs.size()));

        for (auto node : clique)
        {
            if (vNbrs.find(node) != vNbrs.end())
            {
                cand.push_back(node);
            }
        }

        InsertNode(cand, v);
        vector<int> C = clique;
        InsertNode(C, v);

        if (cand == C)
        {
            InsertNode(clique, v);
            cliquesOfNode[v].insert(cliqueIndex);
            resultCliques.insert(cliqueIndex);
        }
        else
        {
            if (MaxEval(cand, cand, G))
            {
                candCliques.insert(move(cand));
            }
        }
    }

    // 添加新的候选cliques
    for (auto &c : candCliques)
    {
        vector<int> clique = c;
        AddClique(maxCliques, cliquesOfNode,
                  notEmptyIndex, emptyIndex, resultCliques, clique);
    }

    // 获取u的邻居集合
    set<int> uNbrs;
    for (const auto &nbr_entry : G[u])
    {
        uNbrs.insert(nbr_entry.first);
    }

    // 处理v所在的cliques
    for (auto cliqueIndex : cliquesV)
    {
        auto &clique = maxCliques[cliqueIndex];
        vector<int> cand;
        cand.reserve(min(clique.size(), uNbrs.size()));

        for (auto node : clique)
        {
            if (uNbrs.find(node) != uNbrs.end())
            {
                cand.push_back(node);
            }
        }

        InsertNode(cand, u);
        vector<int> C = clique;
        InsertNode(C, u);

        if (cand == C)
        {
            DelClique(maxCliques, cliquesOfNode,
                      notEmptyIndex, emptyIndex, resultCliques, cliqueIndex);
        }
    }
}
std::unordered_map<int, std::vector<int>> &GetSubGraph1(
    const std::unordered_map<int, std::unordered_map<int, int>> &graph,
    const std::vector<int> &nodes,
    std::unordered_map<int, std::vector<int>> &subgraph)
{
    if (graph.empty() || nodes.empty())
    {
        return subgraph;
    }

    int nodeCnt = nodes.size();
    int maxNId = nodes.back();
    std::unordered_set<int> nodeSet;
    nodeSet.reserve(nodeCnt);

    for (auto &n : nodes)
    {
        nodeSet.insert(n);
    }

    subgraph.clear();
    // subgraph.reserve(nodeCnt);

    for (auto &n : nodes)
    {
        if (!graph.count(n))
        {
            continue; // 跳过不在图中的节点
        }

        const auto &neighbors = graph.at(n);
        std::vector<int> nbrs;
        nbrs.reserve(neighbors.size());

        for (const auto &[nbr, _] : neighbors)
        {
            if (nbr > maxNId)
            {
                break;
            }
            if (nodeSet.find(nbr) != nodeSet.end())
            {
                nbrs.push_back(nbr);
            }
        }

        subgraph[n] = std::move(nbrs);
    }

    return subgraph;
}
// bool IsIn(const std::vector<int>& smallSet, const std::vector<int>& bigSet)
// {
// 	int first1 = 0, last1 = bigSet.size(), first2 = 0, last2 = smallSet.size();
// 	while (first1 != last1 && first2 != last2)
// 	{
// 		int uNbr = bigSet[first1], vNbr = smallSet[first2];
// 		if (uNbr < vNbr)
// 		{
// 			first1++;
// 		}
// 		else if (uNbr > vNbr)
// 		{
// 			return false;
// 		}
// 		else
// 		{
// 			first1++;
// 			first2++;
// 		}
// 	}
// 	return first2 == last2;
// }

// int SearchNode(const std::vector<int>& clique, int node)
// {
// 	int size = clique.size();
// 	if (size == 0)
// 	{
// 		return -1;
// 	}
// 	return BinSearchNode(clique, node, 0, size - 1);
// }
bool DeleteNode(std::vector<int> &clique, int node)
{
    auto it = find(clique.begin(), clique.end(), node);
    if (it == clique.end())
    {
        return false;
    }
    clique.erase(it);
    // clique.erase(clique.begin() + index);
    return true;
}
vector<vector<int>> AddEdges_NIEMCH(unordered_map<int, unordered_map<int, int>> &G,
                                    vector<vector<int>> &maxCliques,
                                    map<int, set<int>> &cliquesOfNode,
                                    set<int> &notEmptyIndex, set<int> &emptyIndex,
                                    set<int> &resultCliques,
                                    unordered_map<int, unordered_map<int, int>> &edgesToInsert,
                                    chrono::high_resolution_clock::time_point startTime)
{
    vector<pair<int, int>> deg;
    GetNonAscDeg(edgesToInsert, deg);

    // set<int> resultCliques;
    vector<vector<int>> newCliques;

    while (!edgesToInsert.empty())
    {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - startTime;
        if (duration.count() > max_time)
        {
            return vector<vector<int>>();
        }
        auto &deg_u = deg.back();
        deg.pop_back();
        int u = deg_u.second;
        auto &uNeighbors = edgesToInsert[u];
        int uDeg = uNeighbors.size();

        if (uDeg <= 0)
        {
            for (auto &nbr_entry : uNeighbors)
            {
                int nbr = nbr_entry.first;
                AddEdge(G,
                        maxCliques,
                        cliquesOfNode,
                        notEmptyIndex, emptyIndex,
                        resultCliques, u, nbr);
            }
        }
        else
        {
            set<int> N;          // N(u, G)+N(u, G1)
            vector<int> newNbrs; // G1
            vector<int> N_2;     //// N_2 is the actual N(u, G1)
            // N_2.reserve(uDeg);
            // newNbrs.reserve(uDeg);

            for (auto &nbr_entry : uNeighbors)
            {
                int nbr = nbr_entry.first;
                if (G.find(nbr) == G.end())
                {
                    newNbrs.push_back(nbr); // G1
                }
                else
                {
                    N.insert(nbr); // G
                }
            }

            unordered_map<int, vector<int>> nbrGraph;
            vector<int> nbrNodes;

            if (G.find(u) != G.end())
            {
                set<int> N_1; // N(u, G)
                auto &uGNeighbors = G[u];
                // N_1.reserve(uGNeighbors.size());

                for (auto &nbr_entry : uGNeighbors)
                {
                    int nbr = nbr_entry.first;
                    N_1.insert(nbr);
                    N.insert(nbr);
                }

                nbrNodes.reserve(N.size());
                for (auto nbr : N)
                {
                    nbrNodes.push_back(nbr);
                    if (N_1.find(nbr) == N_1.end())
                    {
                        N_2.push_back(nbr);
                    }
                }

                sort(nbrNodes.begin(), nbrNodes.end());
                sort(N_2.begin(), N_2.end());
                nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph); // G[N]

                auto &uCliques = cliquesOfNode[u];
                for (auto it = uCliques.begin(); it != uCliques.end();)
                {
                    int cliqueIndex = *it;
                    auto &clique = maxCliques[cliqueIndex];
                    DeleteNode(clique, u);

                    if (!MaxEval(clique, N_1, nbrGraph))
                    {
                        DelClique(maxCliques,
                                  cliquesOfNode,
                                  notEmptyIndex, emptyIndex,
                                  resultCliques, cliqueIndex);
                        it = uCliques.erase(it);
                    }
                    else
                    {
                        InsertNode(clique, u);
                        ++it;
                    }
                }

                for (int node : N_2)
                {
                    G[u][node] = -1; // 设置边信息为-1
                    G[node][u] = -1;
                }
            }
            else
            {
                // nbrNodes.reserve(N.size());
                for (auto nbr : N)
                {
                    nbrNodes.push_back(nbr);
                }
                sort(nbrNodes.begin(), nbrNodes.end());
                nbrGraph = GetSubGraph1(G, nbrNodes, nbrGraph);

                // 添加新节点u到G中，设置边信息为-1
                for (int nbr : nbrNodes)
                {
                    G[u][nbr] = -1;
                    G[nbr][u] = -1;
                }
                N_2 = std::move(nbrNodes);
            }

            set<int> cliquesToIntersect;
            for (auto n : N_2)
            {
                auto &cliques = cliquesOfNode[n];
                for (auto cliqueIndex : cliques)
                {
                    cliquesToIntersect.insert(cliqueIndex);
                }
            }

            set<vector<int>> candCliques;

            for (auto cliqueIndex : cliquesToIntersect)
            {
                auto &clique = maxCliques[cliqueIndex];
                vector<int> intersections;
                // intersections.reserve(min(clique.size(), N.size()) + 1);

                for (auto node : clique)
                {
                    if (N.find(node) != N.end())
                    {
                        intersections.push_back(node);
                    }
                }

                if (intersections.size() == clique.size())
                {
                    InsertNode(clique, u);
                    cliquesOfNode[u].insert(cliqueIndex);
                    resultCliques.insert(cliqueIndex);
                }
                else
                {
                    candCliques.insert(std::move(intersections));
                }
            }

            for (auto &c : candCliques)
            {
                vector<int> clique = c;
                if (MaxEval(clique, clique, nbrGraph))
                {
                    InsertNode(clique, u);
                    AddClique(maxCliques, cliquesOfNode,
                              notEmptyIndex, emptyIndex, resultCliques, clique);
                }
            }

            for (auto nbr : newNbrs)
            {
                AddEdge(G,
                        maxCliques,
                        cliquesOfNode,
                        notEmptyIndex, emptyIndex,
                        resultCliques, u, nbr);
            }
        }

        edgesToInsert.erase(u);
    }

    // newCliques.reserve(resultCliques.size());
    for (auto cliqueIndex : resultCliques)
    {
        newCliques.push_back(maxCliques[cliqueIndex]);
    }

    resultCliques.clear();
    return move(newCliques);
}
std::unordered_map<int, std::unordered_map<int, int>> GetSubGraph(
    const std::unordered_map<int, std::unordered_map<int, int>> &graph,
    const std::vector<int> &nodes)
{
    std::unordered_map<int, std::unordered_map<int, int>> subgraph;
    std::unordered_set<int> nodeSet(nodes.begin(), nodes.end());

    for (int n : nodes)
    {
        if (!graph.count(n))
            continue;
        subgraph[n]; // init empty
        for (const auto &neighborEntry : graph.at(n))
        {
            int neighbor = neighborEntry.first;
            if (nodeSet.count(neighbor))
            {
                subgraph[n][neighbor] = -1;
                subgraph[neighbor][n] = -1; // undirected
            }
        }
    }
    return move(subgraph);
}
void RemoveReplicas(const std::vector<std::vector<int>> &cliques,
                    std::vector<int> &candIndexes,
                    int candBegin, int candEnd,
                    int level, int size,
                    std::vector<int> &cliqueIndexes)
{
    if (level >= size)
    {
        for (int i = candBegin + 1; i < candEnd + 1; i++)
        {
            int index = candIndexes[i];
            cliqueIndexes[index] = -1;
        }
        return;
    }

    std::vector<int> indexes;
    indexes.reserve(candEnd - candBegin + 1);
    for (int i = candBegin; i < candEnd + 1; i++)
    {
        indexes.push_back(candIndexes[i]);
    }

    std::sort(indexes.begin(), indexes.end(), [&cliques, &level](const int &i1, const int &i2)
              { return cliques[i1][level] < cliques[i2][level]; });

    int left = 0, indexCnt = indexes.size();
    for (int i = 0; i < indexCnt; i++)
    {
        if (cliques[indexes[left]][level] == cliques[indexes[i]][level])
        {
            continue;
        }
        else
        {
            if (i - left >= 2)
            {
                RemoveReplicas(cliques, indexes, left, i - 1, level + 1, size, cliqueIndexes);
            }
            left = i;
        }
    }

    if (left < indexCnt - 1)
    {
        RemoveReplicas(cliques, indexes, left, indexCnt - 1, level + 1, size, cliqueIndexes);
    }
}
unordered_map<int, unordered_map<int, int>> GetMaxKTruss1(
    unordered_map<int, unordered_map<int, int>> &graph,
    int k,
    unordered_map<int, unordered_map<int, int>> &weightEqPairs)
{
    // 检查边数条件
    int edgeCount = 0;
    for (const auto &u_entry : graph)
    {
        edgeCount += u_entry.second.size();
    }
    edgeCount /= 2; // 每条边被计数两次

    int weightEdgeCount = 0;
    for (const auto &u_entry : weightEqPairs)
    {
        weightEdgeCount += u_entry.second.size();
    }
    weightEdgeCount /= 2;

    if (edgeCount < k * (k - 1) / 2 || weightEdgeCount == 0)
    {
        return {};
    }

    // 复制图
    // auto graph = G;
    // map<Edge, int> sup;
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
    // vector<int> cv = {0, 1, 3, 6, 17, 18, 19, 42, 59, 71, 81, 96, 177, 219, 304};
    // vector<int> cv = {3, 6, 17, 18 ,19, 59, 71, 81, 177, 203 ,219, 304};
    // int f = true;
    // for(int u : cv){
    //     if(!graph.count(u)){
    //         f = false;
    //         break;
    //     }
    // }
    // if (f)
    // {
    //     int a = 1;
    //     printGraphAsEdgeList(graph);
    // }
    // 处理边集
    size_t i = 0;
    for (; i < supUqEdges.size(); i++)
    {
        if (supQEdgeCnt == 0 || weightEqPairs.empty())
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
        // 更新weightEqPairs
        if (graph[u].empty() && weightEqPairs.count(u))
        {
            for (auto &v_entry : weightEqPairs[u])
            {
                int v = v_entry.first;
                if (weightEqPairs[v].count(u))
                {
                    weightEqPairs[v].erase(u);
                }
            }
            weightEqPairs.erase(u);
        }
        if (graph[v].empty() && weightEqPairs.count(v))
        {
            for (auto &u_entry : weightEqPairs[v])
            {
                int u = u_entry.first;
                if (weightEqPairs[u].count(v))
                {
                    weightEqPairs[u].erase(v);
                }
            }
            weightEqPairs.erase(u);
        }
    }
    // if(i!=supUqEdges.size()){
    //     int a = 1;
    // }
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
    if (weightEqPairs.empty())
    {
        return {};
    }

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
        weightEqPairs.erase(node);
    }

    return move(graph);
}
unordered_map<int, unordered_map<int, int>> GetMaxKTruss2(
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
// Weight cal_score2(Weight u_score, Weight v_score, Weight u_v_sim)
// {
//     return move((u_score + v_score) * u_v_sim);
// }
bool TestTrussWeightEqual(int method, vector<int> &clique, Weight weight, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    int nodeCnt = clique.size();
    for (int i = 0; i < nodeCnt; i++)
    {
        int u = clique[i];
        Weight u_score = seq2att[u].second.first;
        vector<int> u_att = seq2att[u].second.second;
        for (int j = i + 1; j < nodeCnt; j++)
        {
            int v = clique[j];
            Weight v_score = seq2att[v].second.first;
            vector<int> v_att = seq2att[v].second.second;
            vector<int> com_att;
            set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
            Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
            Weight w = cal_score(u_score, v_score, u_v_sim);
            // Weight w;
            // if (method == 0 || method == 1)
            // {
            //     w = cal_score(u_score, v_score, u_v_sim);
            // }
            // else if (method == 2)
            // {
            //     w = cal_score2(u_score, v_score, u_v_sim);
            // }
            if (w == weight)
            {
                return true;
            }
        }
    }
    return false;
}
void GetResults(int method, int k, unordered_map<int, unordered_map<int, int>> &T,
                Weight weight, vector<vector<int>> &newCliques,
                const unordered_map<int, unordered_map<int, int>> &updateEdgeGraph,
                map<Weight, vector<vector<int>>, greater<Weight>> &results,
                unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
{
    if (weight == Weight(5, 6392))
    {
        int a = 1;
    }
    vector<vector<int>> cliques;
    vector<int> candIndexes, indexes;
    int newCliqueCnt = newCliques.size();
    candIndexes.reserve(newCliqueCnt);
    cliques.reserve(newCliqueCnt);
    indexes.reserve(newCliqueCnt);

    unordered_map<int, vector<int>> cliquesOfSize;

    for (auto &newClique : newCliques)
    {
        // 获取子图
        auto cliqueGraph = GetSubGraph(T, newClique);

        auto weightEqPairs = GetSubGraph(updateEdgeGraph, newClique);

        // 计算最大k-truss
        cliqueGraph = GetMaxKTruss1(cliqueGraph, k, weightEqPairs);
        if (cliqueGraph.empty())
        {
            continue;
        }

        // 获取最大连接的k-truss
        // auto TCT = GetMaxTCKTruss2(cliqueGraph);
        vector<vector<int>> TCT = FindConnectedComponents(cliqueGraph);

        for (auto &clique : TCT)
        {
            if (TestTrussWeightEqual(method, clique, weight, seq2att))
            {
                if (clique.size() == newClique.size())
                {
                    // 已经是最大的
                    results[weight].push_back(clique);
                    cliques.push_back(clique);
                }
                else
                {
                    // line 6
                    int newIndex = cliques.size();
                    cliques.push_back(clique);
                    candIndexes.push_back(newIndex);
                    cliquesOfSize[clique.size()].push_back(newIndex);
                }
                indexes.push_back(indexes.size());
            }
        }
    }

    // 移除相同的连接k-truss
    for (auto &size_indexes : cliquesOfSize)
    {
        int size = size_indexes.first;
        auto &indexes_ = size_indexes.second;
        if (indexes_.size() > 1)
        {
            RemoveReplicas(cliques, indexes_, 0, indexes_.size() - 1, 0, size, indexes);
        }
    }

    vector<pair<int, int>> sortedIndexes;
    sortedIndexes.reserve(indexes.size());
    for (auto index : indexes)
    {
        if (index == -1)
            continue;
        auto &clique = cliques[index];
        sortedIndexes.emplace_back(clique.size(), index);
    }
    sort(sortedIndexes.begin(), sortedIndexes.end(),
         [](const pair<int, int> &a, const pair<int, int> &b)
         {
             return a.first > b.first; // 降序排序
         });

    for (auto index : candIndexes)
    {
        if (indexes[index] == -1)
            continue;

        auto &clique = cliques[index];
        int comNbrCnt = clique.size();
        bool isMaximal = true;

        for (int m = 0; m < sortedIndexes.size(); m++)
        {
            int size = sortedIndexes[m].first;
            if (size <= comNbrCnt)
                break;

            auto bigIndex = sortedIndexes[m].second;
            auto &bigClique = cliques[bigIndex];

            sort(bigClique.begin(), bigClique.end());
            sort(clique.begin(), clique.end());
            if (includes(bigClique.begin(), bigClique.end(),
                         clique.begin(), clique.end()))
            {
                isMaximal = false;
                break;
            }
        }

        if (isMaximal)
        {
            results[weight].push_back(clique);
        }
    }
}
// void QueryByPairs_IEMC(int k, unordered_map<int, unordered_map<int, int>> &T,
//                        unordered_map<int, unordered_map<int, int>> &G,
//                        vector<vector<int>> &maxCliques,
//                        map<int, set<int>> &cliquesOfNode,
//                        set<int> &notEmptyIndex, set<int> &emptyIndex,
//                        set<int> &resultCliques,
//                        const vector<Edge> &pairs,
//                        const vector<pair<Weight, int>> &weightIndex,
//                        map<Weight, vector<vector<int>>, greater<Weight>> &results,
//                        unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att)
// {
//     int pairLIndex = 0;
//     unordered_map<int, unordered_map<int, int>> updateEdgeGraph;

//     for (auto &weight_cnt : weightIndex)
//     {
//         Weight weight = weight_cnt.first;
//         // cout  << "weight: " << weight.toString() << endl;
//         if (weight == Weight(1, 3145))
//         {
//             int a = 1;
//         }
//         int cnt = weight_cnt.second;

//         // 添加边到更新图中（边值初始化为-1）
//         for (int i = pairLIndex; i < pairLIndex + cnt; i++)
//         {
//             int u = pairs[i].first;
//             int v = pairs[i].second;
//             updateEdgeGraph[u][v] = -1;
//             updateEdgeGraph[v][u] = -1;
//         }

//         // 创建图的副本
//         auto updateEdgeGraph1 = updateEdgeGraph;
//         vector<vector<int>> newCliques;

//         // newCliques = dynamicClique->AddEdges_NIEMCH(updateEdgeGraph);
//         newCliques = AddEdges_NIEMCH(G, maxCliques, cliquesOfNode, notEmptyIndex, emptyIndex,
//                                      resultCliques, updateEdgeGraph, start_time);

//         // 获取结果
//         GetResults(k, T, weight, newCliques, updateEdgeGraph1, results, seq2att);
//         pairLIndex += cnt;

//         // 清空图
//         updateEdgeGraph.clear();
//         updateEdgeGraph1.clear();
//     }
// }
// void printGraphAsEdgeList(const std::unordered_map<int, std::unordered_map<int, int>> &graph)
// {
//     std::cout << "edges = [\n";
//     for (auto &u_v : graph)
//     {
//         int src = u_v.first;
//         for (const auto &[dst, weight] : graph.at(src))
//         {
//             std::cout << "    (" << src << ", " << dst << ", " << weight << "),\n";
//         }
//     }
//     std::cout << "]\n";
// }
void record_results(vector<int> cand, Weight &score, vector<pair<Weight, vector<int>>> &result)
{

    sort(cand.begin(), cand.end());

    if (result.size() == 0)
    {
        result.push_back(make_pair(score, cand));
    }
    else
    {
        for (int i = 0; i < result.size(); i++)
        {
            vector<int> &members = result[i].second;
            vector<int> com;
            set_intersection(cand.begin(), cand.end(), members.begin(), members.end(), back_inserter(com));
            if (com.size() == members.size())
            {                   //|com|<=|result[i]|
                members = cand; // cand>>result[i]
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
                        result.push_back(make_pair(score, cand));
                        // cout << "cand score: " << score << endl;
                    }
                }
            }
        }
    }
}
// void NaiveEnum(Weight &score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
//                unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight sthd1, Weight sthd2, vector<pair<Weight, vector<int>>> &result)
// {
//     if (C.empty())
//     {
//         if (M.empty())
//             return;
//         Weight ns = NSize(M, seq2att);
//         if (ns >= sthd1 && ns <= sthd2 && isCntKTruss(M, C, kct, k))
//         {
//             // record
//             // result.push_back(make_pair(score, M));
//             vector<int> cand;
//             for (const auto &it : M)
//             {
//                 cand.push_back(it);
//             }
//             sort(cand.begin(), cand.end());
//             max_compare(cand, score, result);
//         }

//         return;
//     }

//     // 选择 C 中的一个顶点 u
//     int u = C.back();
//     C.pop_back();

//     // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
//     M.push_back(u);

//     NaiveEnum(score, M, C, kct, seq2att, k, sthd1, sthd2, result);

//     M.pop_back(); // 恢复 M

//     unordered_map<int, int> u_neighbors = kct[u];
//     // local update support
//     vector<Edge> changed_edges;
//     set<int> u_neighbors_set;
//     for (auto &v : u_neighbors)
//     {
//         u_neighbors_set.insert(v.first);
//     }
//     for (auto &v : u_neighbors)
//     {
//         set<int> v_neighbors_set;
//         for (auto &it : kct[v.first])
//         {
//             v_neighbors_set.insert(it.first);
//         }
//         set<int> com_neighbors;
//         set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
//         for (int w : com_neighbors)
//         {
//             if (find(changed_edges.begin(), changed_edges.end(), Edge(v.first, w)) == changed_edges.end())
//             {
//                 changed_edges.push_back(Edge(v.first, w));
//                 changed_edges.push_back(Edge(w, v.first));
//                 kct[w][v.first]--;
//                 kct[v.first][w]--;
//             }
//         }
//     }
//     // delete u from kct
//     kct.erase(u);
//     for (const auto &it : u_neighbors)
//     {
//         if (kct.count(it.first))
//         {
//             kct[it.first].erase(u);
//         }
//     }

//     NaiveEnum(score, M, C, kct, seq2att, k, sthd1, sthd2, result);

//     // 恢复 C
//     C.push_back(u);
//     kct.emplace(u, u_neighbors);
//     for (const auto &it : u_neighbors)
//     {
//         if (kct.count(it.first))
//         {
//             kct[it.first].emplace(u, it.second);
//         }
//     }
//     // recover support
//     for (auto &e : changed_edges)
//     {
//         kct[e.first][e.second]++;
//     }
// }
// void NaiveEnum_Truss(Weight &score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
//                      unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight sthd1, Weight sthd2, vector<pair<Weight, vector<int>>> &result)
// {
//     if (C.empty())
//     {
//         // check connectivity
//         unordered_set<int> visited;
//         int startNode = M[0];
//         DFS(kct, startNode, visited);
//         if (visited.size() != M.size())
//         {
//             return; // 图不连通
//         }
//         Weight ns = NSize(M, seq2att);
//         if (ns >= sthd1 && ns <= sthd2)
//         {
//             // record
//             // result.push_back(make_pair(score, M));
//             vector<int> cand;
//             for (const auto &it : M)
//             {
//                 cand.push_back(it);
//             }
//             sort(cand.begin(), cand.end());
//             max_compare(cand, score, result);
//         }

//         return;
//     }

//     // 选择 C 中的一个顶点 u
//     int u = C.back();
//     C.pop_back();

//     // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
//     M.push_back(u);

//     NaiveEnum_Truss(score, M, C, kct, seq2att, k, sthd1, sthd2, result);

//     M.pop_back(); // 恢复 M

//     vector<VertexInfo> changed_edges;
//     set<int> del_vertex;
//     changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
//     // del_vertex.insert(u);

//     set<int> M_temp(M.begin(), M.end());
//     set<int> com;
//     set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
//     if (com.size() == 0)
//     {
//         for (int v : del_vertex)
//         {
//             auto it = std::find(C.begin(), C.end(), v);
//             if (it != C.end())
//             {
//                 // 如果找到，则删除
//                 C.erase(it);
//             }
//         }
//         if (M.size() + C.size() >= k)
//         {
//             NaiveEnum_Truss(score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//         }
//         for (int v : del_vertex)
//             C.push_back(v);
//     }

//     // 恢复 C
//     C.push_back(u);
//     for (int v : del_vertex)
//     {
//         kct.emplace(v, unordered_map<int, int>());
//     }
//     if (!changed_edges.empty())
//         recoverTruss(changed_edges.back(), kct);
// }
// void NaiveEnum_Truss(Weight &score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
//                      unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight &opt, Weight sthd1, Weight sthd2, vector<pair<Weight, vector<int>>> &result)
// {
//     if (C.empty())
//     {
//         // check connectivity
//         unordered_set<int> visited;
//         int startNode = M[0];
//         DFS(kct, startNode, visited);
//         if (visited.size() != M.size())
//         {
//             return; // 图不连通
//         }
//         Weight ns = NSize(M, seq2att);
//         if (ns >= sthd1 && ns <= sthd2)
//         {
//             // record
//             // result.push_back(make_pair(score, M));
//             vector<int> cand;
//             for (const auto &it : M)
//             {
//                 cand.push_back(it);
//             }
//             sort(cand.begin(), cand.end());
//             max_compare(cand, score, opt, result);
//         }

//         return;
//     }

//     // 选择 C 中的一个顶点 u
//     int u = C.back();
//     C.pop_back();

//     // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
//     M.push_back(u);

//     // for (auto &it : M)
//     // {
//     //     cout << it << " ";
//     // }
//     // cout << ", ";
//     // for (auto &it : C)
//     // {
//     //     cout << it << " ";
//     // }
//     // cout << endl;

//     NaiveEnum_Truss(score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
//     M.pop_back(); // 恢复 M
//     // MaintainTruss_DeleteVertex();

//     // delete vertex
//     unordered_map<int, int> u_neighbors = kct[u];
//     // local update support
//     vector<Edge> changed_edges;
//     vector<pair<Edge, int>> sup_less_edges;
//     set<int> u_neighbors_set;
//     for (auto &v : u_neighbors)
//     {
//         u_neighbors_set.insert(v.first);
//     }
//     for (auto &v : u_neighbors)
//     {
//         set<int> v_neighbors_set;
//         for (auto &it : kct[v.first])
//         {
//             v_neighbors_set.insert(it.first);
//         }
//         set<int> com_neighbors;
//         set_intersection(u_neighbors_set.begin(), u_neighbors_set.end(), v_neighbors_set.begin(), v_neighbors_set.end(), inserter(com_neighbors, com_neighbors.begin()));
//         for (int w : com_neighbors)
//         {
//             if (find(changed_edges.begin(), changed_edges.end(), Edge(v.first, w)) == changed_edges.end())
//             {
//                 changed_edges.push_back(Edge(v.first, w));
//                 changed_edges.push_back(Edge(w, v.first));
//                 if (--kct[w][v.first] < (k - 2))
//                 {
//                     sup_less_edges.push_back(make_pair(Edge(w, v.first), 0));
//                 }
//                 kct[v.first][w]--;
//             }
//         }
//     }
//     // delete u from kct
//     kct.erase(u);
//     for (const auto &it : u_neighbors)
//     {
//         if (kct.count(it.first))
//         {
//             kct[it.first].erase(u);
//         }
//     }

//     // delete the edges whose sup < k-2
//     set<int> del_vertex;
//     for (int i = 0; i < sup_less_edges.size(); i++)
//     {
//         Edge e = sup_less_edges[i].first;
//         int u = e.first;
//         int v = e.second;
//         sup_less_edges[i].second = kct[u][v];
//         kct[u].erase(v);
//         kct[v].erase(u);
//         // delete isolated vertices
//         if (kct[u].empty())
//             del_vertex.insert(u);
//         if (kct[v].empty())
//             del_vertex.insert(v);
//         for (auto nei : kct[u])
//         {
//             int w = nei.first;
//             if (kct[v].count(w))
//             {
//                 if (kct[u][w] == (k - 2))
//                 {
//                     sup_less_edges.push_back(make_pair(Edge(u, w), 0));
//                 }
//                 --kct[u][w];
//                 --kct[w][u];
//                 if (kct[v][w] == (k - 2))
//                 {
//                     sup_less_edges.push_back(make_pair(Edge(v, w), 0));
//                 }
//                 --kct[v][w];
//                 --kct[w][v];
//             }
//         }
//     }
//     set<int> M_temp(M.begin(), M.end());
//     set<int> com;
//     set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
//     if (com.size() == 0)
//     {
//         for (int v : del_vertex)
//         {
//             auto it = std::find(C.begin(), C.end(), v);
//             if (it != C.end())
//             {
//                 // 如果找到，则删除
//                 C.erase(it);
//             }
//         }
//         if (M.size() + C.size() >= k)
//         {
//             NaiveEnum_Truss(score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
//         }
//         for (int v : del_vertex)
//             C.push_back(v);
//     }

//     // 恢复 C
//     C.push_back(u);
//     while (!sup_less_edges.empty())
//     {
//         auto &it = sup_less_edges.back();

//         Edge e = it.first;
//         int u = e.first;
//         int v = e.second;
//         for (auto nei : kct[u])
//         {
//             int w = nei.first;
//             if (kct[v].count(w))
//             {
//                 ++kct[u][w];
//                 ++kct[w][u];
//                 ++kct[v][w];
//                 ++kct[w][v];
//             }
//         }
//         kct[u].emplace(v, it.second);
//         kct[v].emplace(u, it.second);
//         sup_less_edges.pop_back();
//     }
//     kct.emplace(u, u_neighbors);
//     for (const auto &it : u_neighbors)
//     {
//         if (kct.count(it.first))
//         {
//             kct[it.first].emplace(u, it.second);
//         }
//     }

//     // recover support
//     for (auto &e : changed_edges)
//     {
//         kct[e.first][e.second]++;
//     }
// }
// void NaiveEnum_CntTruss(Weight &score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
//                         unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight sthd1, Weight sthd2, vector<pair<Weight, vector<int>>> &result)
// {
//     if (C.empty())
//     {
//         // check connectivity
//         unordered_set<int> visited;
//         int startNode = M[0];
//         DFS(kct, startNode, visited);
//         if (visited.size() != M.size())
//         {
//             return; // 图不连通
//         }
//         Weight ns = NSize(M, seq2att);
//         if (ns >= sthd1 && ns <= sthd2)
//         {
//             // record
//             // result.push_back(make_pair(score, M));
//             vector<int> cand;
//             for (const auto &it : M)
//             {
//                 cand.push_back(it);
//             }
//             sort(cand.begin(), cand.end());
//             max_compare(cand, score, result);
//         }

//         return;
//     }

//     // 选择 C 中的一个顶点 u
//     int u = C.back();
//     C.pop_back();

//     // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
//     M.push_back(u);

//     NaiveEnum_CntTruss(score, M, C, kct, seq2att, k, sthd1, sthd2, result);

//     M.pop_back(); // 恢复 M

//     vector<VertexInfo> changed_edges;
//     set<int> del_vertex;
//     changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
//     // del_vertex.insert(u);

//     set<int> M_temp(M.begin(), M.end());
//     set<int> com;
//     set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
//     if (com.size() == 0)
//     {
//         for (int v : del_vertex)
//         {
//             auto it = std::find(C.begin(), C.end(), v);
//             if (it != C.end())
//             {
//                 // 如果找到，则删除
//                 C.erase(it);
//             }
//         }
//         // if (M.size() + C.size() >= k)
//         // {
//         //     NaiveEnum_CntTruss(score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
//         // }
//         if (M.empty())
//         {
//             vector<vector<int>> components = FindConnectedComponents(kct, C);
//             if (components.size() == 1)
//             {
//                 NaiveEnum_CntTruss(score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//             }
//             else
//             {
//                 for (auto &c : components)
//                 {
//                     if (c.empty())
//                     {
//                         continue;
//                     }
//                     unordered_map<int, unordered_map<int, int>> local_kct;
//                     for (auto v : c)
//                     {
//                         local_kct.emplace(v, kct.at(v));
//                     }
//                     NaiveEnum_CntTruss(score, M, c, local_kct, seq2att, k, sthd1, sthd2, result);
//                 }
//             }
//         }
//         else
//         {

//             unordered_set<int> visited;
//             int startNode = M[0];
//             DFS(kct, startNode, visited);
//             if (visited.size() == M.size() + C.size())
//             {
//                 NaiveEnum_CntTruss(score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//             }
//             else
//             {
//                 vector<int> C_temp;
//                 int M_count = 0;
//                 for (int v : visited)
//                 {
//                     if (count(M.begin(), M.end(), v))
//                         M_count++;
//                     else
//                         C_temp.push_back(v);
//                 }
//                 if (M_count == M.size())
//                 {
//                     unordered_map<int, unordered_map<int, int>> local_kct;
//                     for (auto v : visited)
//                     {
//                         local_kct.emplace(v, kct.at(v));
//                     }
//                     NaiveEnum_CntTruss(score, M, C_temp, local_kct, seq2att, k, sthd1, sthd2, result);
//                 }
//             }
//         }
//         for (int v : del_vertex)
//             C.push_back(v);
//     }

//     // 恢复 C
//     C.push_back(u);
//     for (int v : del_vertex)
//     {
//         kct.emplace(v, unordered_map<int, int>());
//     }
//     if (!changed_edges.empty())
//         recoverTruss(changed_edges.back(), kct);
// }
bool CConNewVertex(vector<int> C, set<int> &new_vertex)
{
    sort(C.begin(), C.end());
    vector<int> com;
    set_intersection(C.begin(), C.end(), new_vertex.begin(), new_vertex.end(), back_inserter(com));
    if (com.empty())
    {
        return false;
    }
    else
    {
        return true;
    }
}
void NaiveEnum_CntTrussNsize(Weight &score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
                             unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight sthd1, Weight sthd2,
                             map<Weight, vector<vector<int>>, greater<Weight>> &result, chrono::high_resolution_clock::time_point startTime)
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = now - startTime;
    if (duration.count() > max_time)
    { // 超过1秒，终止递归
        // std::cerr << "Timeout: Execution exceeded " << max_time << " second." << std::endl;
        result.clear();
        return;
    }
    // cout << "M: ";
    // for (int u : M)
    // {
    //     cout << u << " ";
    // }
    // cout << endl;
    // cout << "C: ";
    // for (int u : C)
    // {
    //     cout << u << " ";
    // }
    // cout << endl;
    // vector<int> m = {177, 59};
    // if (M == m)
    // {
    //     int a = 1;
    //     printGraphAsEdgeList(kct);
    // }
    if (C.empty())
    {
        // // check connectivity
        // unordered_set<int> visited;
        // int startNode = M[0];
        // DFS(kct, startNode, visited);
        // if (visited.size() != M.size())
        // {
        //     return; // 图不连通
        // }
        Weight ns = NSize(M, seq2att);
        if (ns >= sthd1 && ns <= sthd2)
        {
            // record
            // result.push_back(make_pair(score, M));
            vector<int> cand;
            for (const auto &it : M)
            {
                cand.push_back(it);
            }
            sort(cand.begin(), cand.end());
            max_compare(cand, score, result);
        }

        return;
    }

    // 选择 C 中的一个顶点 u
    int u = C.back();
    C.pop_back();

    // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
    M.push_back(u);
    if (CheckNsize(M, C, seq2att, sthd1, sthd2))
    {
        // if(MConNewVert || CConNewVertex(C, new_vertex))
        NaiveEnum_CntTrussNsize(score, M, C, kct, seq2att, k, sthd1, sthd2, result, startTime);
    }

    M.pop_back(); // 恢复 M

    vector<VertexInfo> changed_edges;
    set<int> del_vertex;
    changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
    // del_vertex.insert(u);

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
        // if (M.size() + C.size() >= k)
        // {
        //     NaiveEnum_CntTruss(score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
        // }
        if (M.empty())
        {
            vector<vector<int>> components = FindConnectedComponents(kct, C);
            if (components.size() == 1)
            {
                NaiveEnum_CntTrussNsize(score, M, C, kct, seq2att, k, sthd1, sthd2, result, startTime);
            }
            else
            {
                for (auto &c : components)
                {
                    if (c.empty())
                    {
                        continue;
                    }
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : c)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_CntTrussNsize(score, M, c, local_kct, seq2att, k, sthd1, sthd2, result, startTime);
                }
            }
        }
        else
        {

            unordered_set<int> visited;
            int startNode = M[0];
            DFS(kct, startNode, visited);
            if (visited.size() == M.size() + C.size() && CheckNsize(M, C, seq2att, sthd1, sthd2))
            {
                NaiveEnum_CntTrussNsize(score, M, C, kct, seq2att, k, sthd1, sthd2, result, startTime);
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
                if (M_count == M.size() && CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
                {
                    unordered_map<int, unordered_map<int, int>> local_kct;
                    for (auto v : visited)
                    {
                        local_kct.emplace(v, kct.at(v));
                    }
                    NaiveEnum_CntTrussNsize(score, M, C_temp, local_kct, seq2att, k, sthd1, sthd2, result, startTime);
                }
            }
        }
        for (int v : del_vertex)
            C.push_back(v);
    }

    // 恢复 C
    C.push_back(u);
    for (int v : del_vertex)
    {
        kct.emplace(v, unordered_map<int, int>());
    }
    if (!changed_edges.empty())
        recoverTruss(changed_edges.back(), kct);
}
// void NaiveEnum_CntTrussNsize(set<int> &new_vertex, Weight &score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
//                              unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight sthd1, Weight sthd2, map<Weight, vector<vector<int>>, greater<Weight>> &result)
// {
//     cout << "new_vertex: " << endl;

//     for (int v : new_vertex)
//     {
//         cout << v << " ";
//     }
//     cout << endl;
//     cout << "M: " << endl;
//     for (int v : M)
//     {
//         cout << v << " ";
//     }
//     cout << endl;
//     cout << "C: " << endl;
//     for (int v : C)
//     {
//         cout << v << " ";
//     }
//     cout << endl;
//     if (C.empty())
//     {
//         // // check connectivity
//         // unordered_set<int> visited;
//         // int startNode = M[0];
//         // DFS(kct, startNode, visited);
//         // if (visited.size() != M.size())
//         // {
//         //     return; // 图不连通
//         // }
//         Weight ns = NSize(M, seq2att);
//         if (ns >= sthd1 && ns <= sthd2)
//         {
//             // record
//             // result.push_back(make_pair(score, M));
//             vector<int> cand;
//             for (const auto &it : M)
//             {
//                 cand.push_back(it);
//             }
//             sort(cand.begin(), cand.end());
//             max_compare(cand, score, result);
//         }

//         return;
//     }

//     // vector<int> C_= {23, 266, 378, 548};
//     // if( C == C_){
//     //     int a = 1;
//     // }

//     // 选择 C 中的一个顶点 u
//     int u = C.back();
//     C.pop_back();

//     // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
//     M.push_back(u);
//     if (CheckNsize(M, C, seq2att, sthd1, sthd2))
//     {
//         // if(MConNewVert || CConNewVertex(C, new_vertex))
//         NaiveEnum_CntTrussNsize(new_vertex, score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//     }

//     M.pop_back(); // 恢复 M

//     vector<VertexInfo> changed_edges;
//     set<int> del_vertex;
//     changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
//     // del_vertex.insert(u);

//     set<int> M_temp(M.begin(), M.end());
//     vector<int> del_new_vertex;
//     set<int> com;
//     set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
//     if (com.size() == 0)
//     {
//         if (new_vertex.count(u))
//         {
//             new_vertex.erase(u);
//             del_new_vertex.push_back(u);
//         }
//         for (int v : del_vertex)
//         {
//             auto it = std::find(C.begin(), C.end(), v);
//             if (it != C.end())
//             {
//                 // 如果找到，则删除
//                 C.erase(it);
//             }
//             if (new_vertex.count(v))
//             {
//                 new_vertex.erase(v);
//                 del_new_vertex.push_back(v);
//             }
//         }
//         if (M.empty())
//         {
//             vector<vector<int>> components = FindConnectedComponents(kct, C);
//             if (components.size() == 1)
//             {
//                 if (!new_vertex.empty())
//                 {
//                     NaiveEnum_CntTrussNsize(new_vertex, score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//                 }
//             }
//             else
//             {
//                 for (auto &c : components)
//                 {
//                     if (c.empty())
//                     {
//                         continue;
//                     }
//                     set<int> local_new_vertex;
//                     sort(c.begin(), c.end());
//                     set_intersection(c.begin(), c.end(), new_vertex.begin(), new_vertex.end(), inserter(local_new_vertex, local_new_vertex.begin()));
//                     if (local_new_vertex.empty())
//                         continue;
//                     unordered_map<int, unordered_map<int, int>> local_kct;
//                     for (auto v : c)
//                     {
//                         local_kct.emplace(v, kct.at(v));
//                     }
//                     NaiveEnum_CntTrussNsize(local_new_vertex, score, M, c, local_kct, seq2att, k, sthd1, sthd2, result);
//                 }
//             }
//         }
//         else
//         {

//             unordered_set<int> visited;
//             int startNode = M[0];
//             DFS(kct, startNode, visited);
//             if (visited.size() == M.size() + C.size() && CheckNsize(M, C, seq2att, sthd1, sthd2) && !new_vertex.empty())
//             {
//                 NaiveEnum_CntTrussNsize(new_vertex, score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//             }
//             else
//             {
//                 vector<int> C_temp;
//                 set<int> new_vertex_temp;
//                 int M_count = 0;
//                 for (int v : visited)
//                 {
//                     if (count(M.begin(), M.end(), v))
//                         M_count++;
//                     else
//                     {
//                         C_temp.push_back(v);
//                         if (new_vertex.count(v))
//                         {
//                             new_vertex_temp.insert(v);
//                         }
//                     }
//                 }
//                 if (M_count == M.size() && CheckNsize(M, C_temp, seq2att, sthd1, sthd2) && !new_vertex_temp.empty())
//                 {
//                     unordered_map<int, unordered_map<int, int>> local_kct;
//                     for (auto v : visited)
//                     {
//                         local_kct.emplace(v, kct.at(v));
//                     }
//                     NaiveEnum_CntTrussNsize(new_vertex_temp, score, M, C_temp, local_kct, seq2att, k, sthd1, sthd2, result);
//                 }
//             }
//         }
//         for (int v : del_new_vertex)
//         {
//             new_vertex.insert(v);
//         }
//         for (int v : del_vertex)
//             C.push_back(v);
//     }

//     // 恢复 C
//     C.push_back(u);
//     for (int v : del_vertex)
//     {
//         kct.emplace(v, unordered_map<int, int>());
//     }
//     if (!changed_edges.empty())
//         recoverTruss(changed_edges.back(), kct);
// }
// int VertexScore(unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
//                 unordered_map<int, unordered_map<int, int>> &kct, vector<int> &M,
//                 vector<int> &CC, Weight &thd1, Weight &thd2)
// {
//     int n_M = 0;
//     for (auto it : M)
//     {
//         if (!seq2att[it].first)
//         {
//             n_M++;
//         }
//     }
//     Weight ns;
//     if (M.empty())
//     {
//         ns = WEIGHT_ZERO;
//     }
//     else
//     {
//         ns = Weight(n_M, M.size());
//     }

//     // vector<int> vert_sort;
//     int max_deg_vert = -1;
//     int max_deg = -1;
//     // double max_score = INT16_MAX;
//     // unordered_map<int, double> vert_score;
//     int neg = 0;
//     if (ns < thd1) // select neg
//     {
//         neg = -1;
//     }
//     else if (ns > thd2) // select pos
//     {
//         neg = 1;
//     }
//     for (int u : CC)
//     {
//         if (neg == -1 && seq2att[u].first)
//             continue;
//         else if (neg == 1 && !seq2att[u].first)
//             continue;

//         int deg = kct[u].size();
//         if (deg > max_deg)
//         {
//             max_deg = deg;
//             max_deg_vert = u;
//         }
//     }
//     // if(max_score_vert == -1){
//     //     int a = 1;
//     // }
//     return max_deg_vert;
// }
// void NaiveEnum_CntTrussNsizeVert(Weight &score, vector<int> &M, vector<int> &C, unordered_map<int, unordered_map<int, int>> &kct,
//                                  unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att, const int k, Weight sthd1, Weight sthd2, vector<pair<Weight, vector<int>>> &result)
// {
//     if (C.empty())
//     {
//         // check connectivity
//         unordered_set<int> visited;
//         int startNode = M[0];
//         DFS(kct, startNode, visited);
//         if (visited.size() != M.size())
//         {
//             return; // 图不连通
//         }
//         Weight ns = NSize(M, seq2att);
//         if (ns >= sthd1 && ns <= sthd2)
//         {
//             // record
//             // result.push_back(make_pair(score, M));
//             vector<int> cand;
//             for (const auto &it : M)
//             {
//                 cand.push_back(it);
//             }
//             sort(cand.begin(), cand.end());
//             max_compare(cand, score, result);
//         }

//         return;
//     }

//     // 选择 C 中的一个顶点 u
//     // int u = C.back();
//     // C.pop_back();
//     int u = VertexScore(seq2att, kct, M, C, sthd1, sthd2); // 可以优化，这里每个搜索点都计算一次，但很多分支的kct不变，导致重复计算
//     C.erase(std::find(C.begin(), C.end(), u));

//     // Expand: 调用 NaiveEnum(M ∪ u, C \ u)
//     M.push_back(u);
//     if (CheckNsize(M, C, seq2att, sthd1, sthd2))
//     {
//         NaiveEnum_CntTrussNsizeVert(score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//     }

//     M.pop_back(); // 恢复 M

//     vector<VertexInfo> changed_edges;
//     set<int> del_vertex;
//     changed_edges.push_back(MaintainTruss_DeleteVertex(u, kct, k, del_vertex));
//     // del_vertex.insert(u);

//     set<int> M_temp(M.begin(), M.end());
//     set<int> com;
//     set_intersection(M_temp.begin(), M_temp.end(), del_vertex.begin(), del_vertex.end(), std::inserter(com, com.begin()));
//     if (com.size() == 0)
//     {
//         for (int v : del_vertex)
//         {
//             auto it = std::find(C.begin(), C.end(), v);
//             if (it != C.end())
//             {
//                 // 如果找到，则删除
//                 C.erase(it);
//             }
//         }
//         // if (M.size() + C.size() >= k)
//         // {
//         //     NaiveEnum_CntTruss(score, M, C, kct, seq2att, k, opt, sthd1, sthd2, result);
//         // }
//         if (M.empty())
//         {
//             vector<vector<int>> components = FindConnectedComponents(kct, C);
//             if (components.size() == 1)
//             {
//                 NaiveEnum_CntTrussNsizeVert(score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//             }
//             else
//             {
//                 for (auto &c : components)
//                 {
//                     if (c.empty())
//                     {
//                         continue;
//                     }
//                     unordered_map<int, unordered_map<int, int>> local_kct;
//                     for (auto v : c)
//                     {
//                         local_kct.emplace(v, kct.at(v));
//                     }
//                     NaiveEnum_CntTrussNsizeVert(score, M, c, local_kct, seq2att, k, sthd1, sthd2, result);
//                 }
//             }
//         }
//         else
//         {

//             unordered_set<int> visited;
//             int startNode = M[0];
//             DFS(kct, startNode, visited);
//             if (visited.size() == M.size() + C.size() && CheckNsize(M, C, seq2att, sthd1, sthd2))
//             {
//                 NaiveEnum_CntTrussNsizeVert(score, M, C, kct, seq2att, k, sthd1, sthd2, result);
//             }
//             else
//             {
//                 vector<int> C_temp;
//                 int M_count = 0;
//                 for (int v : visited)
//                 {
//                     if (count(M.begin(), M.end(), v))
//                         M_count++;
//                     else
//                         C_temp.push_back(v);
//                 }
//                 if (M_count == M.size() && CheckNsize(M, C_temp, seq2att, sthd1, sthd2))
//                 {
//                     unordered_map<int, unordered_map<int, int>> local_kct;
//                     for (auto v : visited)
//                     {
//                         local_kct.emplace(v, kct.at(v));
//                     }
//                     NaiveEnum_CntTrussNsizeVert(score, M, C_temp, local_kct, seq2att, k, sthd1, sthd2, result);
//                 }
//             }
//         }
//         for (int v : del_vertex)
//             C.push_back(v);
//     }

//     // 恢复 C
//     C.push_back(u);
//     for (int v : del_vertex)
//     {
//         kct.emplace(v, unordered_map<int, int>());
//     }
//     if (!changed_edges.empty())
//         recoverTruss(changed_edges.back(), kct);
// }
// unordered_map<int, pair<bool, pair<Weight, vector<int>>>> LoadAtt(string dataset, vector<int> posKws)
// {
//     unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att;
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
//         vector<int> comKw;
//         set_intersection(posKws.begin(), posKws.end(), att.begin(), att.end(), back_inserter(comKw));
//         if (comKw.size() > 0)
//         {

//             pair<Weight, vector<int>> kscoreAtt(Weight(comKw.size(), 1), att);
//             id2att.emplace(ver, pair<bool, pair<Weight, vector<int>>>(false, kscoreAtt));
//         }
//     }
//     ain.close();
//     cout << "Loaded attribute successfully!" << endl;
//     return move(id2att);
// }
/*method 0: size + n, 1: not size + n, 2: no size + p*/
map<Weight, vector<vector<int>>, greater<Weight>> Best_First(int method, string dataset, vector<int> posKws, vector<int> negKws, 
    int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    // auto start_time = std::chrono::high_resolution_clock::now();
    map<Weight, vector<vector<int>>, greater<Weight>> size_results;
    // auto start = std::chrono::high_resolution_clock::now();
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att;
    int resultCnt = 0;
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    if (method == 0 || method == 1)
    {
        id2att = LoadAtt(0, AttFileName, posKws, negKws);
    }
    else if (method == 2)
    {
        id2att = LoadAtt(1, AttFileName, posKws, negKws);
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    // find connected ktruss
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    if (kct.empty())
    {
        cout << "No k-truss found!" << endl;
        return size_results;
    }
    else{
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
    for (auto &it1 : kct)
    {
        int u = it1.first;
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
            // Weight u_v;
            // if (method == 0 || method == 1)
            // {
            //     u_v = cal_score(u_score, v_score, u_v_sim);
            // }
            // else if (method == 2)
            // {
            //     u_v = cal_score2(u_score, v_score, u_v_sim);
            // }

            if (!score_edges.count(u_v))
            {
                score_edges.emplace(u_v, vector<pair<int, int>>());
                // score_edges.emplace(u_v, unordered_map<int, set<int>>());
            }
            score_edges[u_v].push_back(pair<int, int>(u, v));
        }
    }
    // cout << "calculate similarity complete. " << endl;
    // map<Weight, vector<vector<int>>, greater<Weight>> results;
    Weight ub;
    if (method == 0 || method == 1)
    {
        ub = Weight(2, 1);
    }
    else if (method == 2)
    {
        ub = Weight(2 * posKws.size(), 1);
    }
    // Weight ub(2, 1);
    map<Weight, vector<Edge>, greater<Weight>> pairsLeqUb;
    unordered_map<int, unordered_map<int, int>> S;
    vector<set<Edge>> nonDecSup; //[1] dont contain the old k-truss edges.
    nonDecSup.resize(2);
    auto it = score_edges.begin();

    unordered_map<int, unordered_map<int, int>> T;
    vector<unordered_map<int, unordered_map<int, int>>> maxCTrusses;
    map<Edge, int> CTrussOfEdge;
    map<int, vector<int>> newPairs;
    queue<int> emptyIndexT;

    unordered_map<int, unordered_map<int, int>> updateEdgeGraph;

    unordered_map<int, unordered_map<int, int>> G;
    vector<vector<int>> maxCliques;                     // maximal cliques of G
    map<int, set<int>> cliquesOfNode;                   // the indices of maximal cliques containing each node in maxCliques
    set<int> notEmptyIndex, emptyIndexC, resultCliques; // the non-empty indices in maxCliques; the empty indices in maxCliques
    int cnt = 0;
    set<int> old_vertex;
    while (it != score_edges.end() || !pairsLeqUb.empty())
    {
        // if (!G.empty())
        // {
        //     printGraphAsEdgeList(G);
        // }

        // for (auto &u_v : G)
        // {
        //     int u = u_v.first;
        //     for (auto &v_s : u_v.second)
        //     {
        //         int v = v_s.first;
        //         if (u == v)
        //         {
        //             int a = 1;
        //         }
        //     }
        // }

        // cout << "cnt: " << cnt << endl;
        // if (cnt == 82)
        // {
        //     int a = 1;
        // }
        cnt++;

        Weight w = it != score_edges.end() ? it->first : WEIGHT_ZERO;

        vector<Edge> P_1;
        vector<Edge> P_2;
        vector<pair<Weight, int>> P_1_WeightIndex;
        if (it != score_edges.end())
        {
            for (auto &e : it->second)
            {
                MaintainSupport_AddEdge(S, nonDecSup, e, k);
            }
            set<Edge> newEdges = GetMaxKTrussInc(T, k, nonDecSup);
            // if (w == Weight(5, 182))
            // {
            //     printGraphAsEdgeList(T);
            //     int a = 1;
            // }
            if (!newEdges.empty())
            {
                AddEdges(newEdges, T,
                         maxCTrusses,
                         CTrussOfEdge, newPairs,
                         emptyIndexT);

                set<int> nodeAttrs;
                // compute the pairwise weight for each pair in P
                for (auto &pair : newPairs)
                {
                    int u = pair.first;
                    Weight u_score = seq2att[u].second.first;
                    vector<int> u_att = seq2att[u].second.second;
                    auto &nodes = pair.second;
                    for (int v : nodes)
                    {
                        Weight v_score = seq2att[v].second.first;
                        vector<int> v_att = seq2att[v].second.second;
                        vector<int> com_att;
                        set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
                        Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
                        Weight weight = cal_score(u_score, v_score, u_v_sim);
                        // Weight weight;
                        // if (method == 0 || method == 1)
                        // {
                        //     weight = cal_score(u_score, v_score, u_v_sim);
                        // }
                        // else if (method == 2)
                        // {
                        //     weight = cal_score2(u_score, v_score, u_v_sim);
                        // }
                        if (weight == WEIGHT_ZERO)
                        {
                            continue;
                        }

                        if (weight < ub)
                        {
                            pairsLeqUb[weight].push_back(make_pair(u, v));
                        }
                        else
                        {
                            P_2.push_back(make_pair(u, v));
                        }
                    }
                    nodeAttrs.clear();
                }
                newPairs.clear();
            }
            it++;
        }
        // cout << "add edges complete." << endl;
        // obtain P_1 and P_0(= P_0 \ P_1)
        for (auto weight_edges = pairsLeqUb.begin(); weight_edges != pairsLeqUb.end();)
        {
            Weight weight = weight_edges->first;
            if (weight >= w && weight < ub)
            {
                P_1.insert(P_1.end(), weight_edges->second.begin(), weight_edges->second.end());
                P_1_WeightIndex.push_back(make_pair(weight, weight_edges->second.size()));
                pairsLeqUb.erase(weight_edges++);
            }
            else
            {
                weight_edges++;
            }
        }

        for (auto &pair : P_2)
        {
            int u = pair.first;
            int v = pair.second;
            if (!updateEdgeGraph.count(u))
            {
                updateEdgeGraph.emplace(u, unordered_map<int, int>());
            }
            updateEdgeGraph[u][v] = 0;
            if (!updateEdgeGraph.count(v))
            {
                updateEdgeGraph.emplace(v, unordered_map<int, int>());
            }
            updateEdgeGraph[v][u] = 0;
        }
        // cout << "begin to niemch" << endl;
        AddEdges_NIEMCH(G, maxCliques, cliquesOfNode, notEmptyIndex, emptyIndexC,
                        resultCliques, updateEdgeGraph, start_time);
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - start_time;
        if (duration.count() > max_time)
        {
            return size_results;
        }
        // cout << "niemch complete" << endl;
        updateEdgeGraph.clear();
        vector<Edge>().swap(P_2);

        // lines 14-18 of Algorithm 6
        // QueryByPairs_IEMC(k, T, G, maxCliques, cliquesOfNode, notEmptyIndex, emptyIndexC,
        //                   resultCliques, P_1, P_1_WeightIndex, results, seq2att);

        int pairLIndex = 0;
        // unordered_map<int, unordered_map<int, int>> updateEdgeGraph;
        bool flag = false;
        for (auto &weight_cnt : P_1_WeightIndex)
        {
            Weight weight = weight_cnt.first;
            // cout << "weight: " << weight.toString() << endl;
            if (weight == Weight(5, 182))
            {
                int a = 1;
            }
            int cnt = weight_cnt.second;

            // 添加边到更新图中（边值初始化为-1）
            for (int i = pairLIndex; i < pairLIndex + cnt; i++)
            {
                int u = P_1[i].first;
                int v = P_1[i].second;
                updateEdgeGraph[u][v] = -1;
                updateEdgeGraph[v][u] = -1;
            }

            // 创建图的副本
            auto updateEdgeGraph1 = updateEdgeGraph;
            vector<vector<int>> newCliques;

            // newCliques = dynamicClique->AddEdges_NIEMCH(updateEdgeGraph);
            newCliques = AddEdges_NIEMCH(G, maxCliques, cliquesOfNode, notEmptyIndex, emptyIndexC,
                                         resultCliques, updateEdgeGraph, start_time);
            // cout << "update niemcn complete" << endl;
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = now - start_time;
            if (duration.count() > max_time)
            {
                return size_results;
            }
            // 获取结果
            map<Weight, vector<vector<int>>, greater<Weight>> results;
            GetResults(method, k, T, weight, newCliques, updateEdgeGraph1, results, seq2att);
            // cout << "get candidates complete" << endl;
            pairLIndex += cnt;

            // vector<set<int>> newVertex;
            // for (auto &w_r : results)
            // {

            //     for (auto &r : w_r.second)
            //     {
            //         set<int> new_vertex;
            //         sort(r.begin(), r.end());
            //         set_difference(r.begin(), r.end(), old_vertex.begin(), old_vertex.end(), inserter(new_vertex, new_vertex.begin()));
            //         old_vertex.insert(new_vertex.begin(), new_vertex.end());
            //         newVertex.push_back(new_vertex);
            //     }
            // }
            if (method == 0)
            {
                // map<Weight, vector<vector<int>>, greater<Weight>> size_results;
                for (auto &res : results)
                {
                    Weight score = res.first;
                    for (int i = 0; i < res.second.size(); i++)
                    {
                        // get subgraph
                        auto &clique = res.second[i];
                        vector<int> m = {3, 25, 168, 177};
                        bool f = true;
                        for (int i : m)
                        {
                            if (!count(clique.begin(), clique.end(), i))
                            {
                                f = false;
                                break;
                            }
                        }
                        if (f)
                        {
                            int a = 1;
                        }
                        unordered_map<int, unordered_map<int, int>> kct = GetSubGraph(T, clique);
                        kct = GetMaxKTruss2(kct, k);
                        unordered_set<int> visited;
                        DFS(kct, clique[0], visited);
                        if (kct.size() != clique.size() || visited.size() != clique.size())
                        {
                            cout << "error" << endl;
                        }
                        vector<int> M;
                        // for (auto &u_v : kct)
                        // {
                        //     int u = u_v.first;
                        //     for (auto &v_s : u_v.second)
                        //     {
                        //         int v = v_s.first;
                        //         if (u < v)
                        //             continue;
                        //         int sup = 0;
                        //         for (auto &w_s : kct[v])
                        //         {
                        //             int w = w_s.first;
                        //             if (u_v.second.count(w))
                        //             {
                        //                 sup++;
                        //             }
                        //         }
                        //         kct[u][v] = sup;
                        //         kct[v][u] = sup;
                        //     }
                        // }
                        // vector<int> m = {3, 6, 17, 18 ,19, 59, 71, 81, 177, 203 ,219, 304};
                        // bool f = true;
                        // for(int i : m){
                        //     if(!count(clique.begin(), clique.end(), i)){
                        //         f = false;
                        //         break;
                        //     }
                        // }
                        // if(f){
                        //     int a = 1;
                        // }
                        // if (count(clique.begin(), clique.end(), 177) && count(clique.begin(), clique.end(), 59) && count(clique.begin(), clique.end(), 1))
                        // {
                        //     int a = 1;
                        //     for(int u  : clique){
                        //         cout << u << ", ";
                        //     }
                        //     cout << endl;
                        //     printGraphAsEdgeList(kct);
                        // }
                        // NaiveEnum(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        // NaiveEnum_Truss(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        // NaiveEnum_CntTruss(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        NaiveEnum_CntTrussNsize(score, M, clique, kct, seq2att, k, sthd1, sthd2, size_results, start_time);
                        // NaiveEnum_CntTrussNsize(newVertex[i], score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        // NaiveEnum_CntTrussNsizeVert(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                    }
                }
                if (!size_results.empty())
                {
                    flag = true;
                    printf_results(dg, size_results, seq2att);
                    break;
                }
            }
            else
            {
                if (!results.empty())
                {

                    flag = true;
                    printf_results(dg, results, seq2att);
                    size_results = results;
                    break;
                }
            }
            results.clear();

            // 清空图
            updateEdgeGraph.clear();
            updateEdgeGraph1.clear();
        }

        ub = w;

        // resultCnt += results.size();
        // if (resultCnt > 10)
        // {
        //     break;
        // }
        if (flag)
        {
            // printf_results(size_results, seq2att);
            break;
        }
        // results.clear();
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;
    // for (auto &res : results)
    // {
    //     cout << "Weight: " << res.first.toString() << ", ";
    //     for (auto &clique : res.second)
    //     {
    //         for (auto &node : clique)
    //         {
    //             cout << node << " ";
    //         }
    //         cout << endl;
    //     }
    // }
    return move(size_results);
}

map<Weight, vector<vector<int>>, greater<Weight>> Best_First_S(int method, string dataset, vector<int> posKws, vector<int> negKws, 
    int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time)
{
    // auto start_time = std::chrono::high_resolution_clock::now();
    map<Weight, vector<vector<int>>, greater<Weight>> size_results;
    // auto start = std::chrono::high_resolution_clock::now();
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att;
    int resultCnt = 0;
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    if (method == 0 || method == 1)
    {
        id2att = LoadAtt(0, AttFileName, posKws, negKws);
    }
    else if (method == 2)
    {
        id2att = LoadAtt(1, AttFileName, posKws, negKws);
    }
    DataGraph dg = LoadGraph(dataset, id2att);
    // find connected ktruss
    TrussDecomposition *td = new TrussDecomposition(dg, k, start_time);
    unordered_map<int, unordered_map<int, int>> &kct = td->kTruss;
    if (kct.empty())
    {
        cout << "No k-truss found!" << endl;
        return size_results;
    }
    else{
        cout << "k-truss found!" << endl;
    }
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> seq2att;
    // map<Weight, vector<pair<int, int>>, CompareKeys> score_edges;
    for (auto &it : kct)
    {
        int idx = dg.seq2id[it.first];
        seq2att[it.first] = id2att[idx];
    }
    id2att.clear();

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
    std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start_time;
    cout << "Time used: " << duration.count() << "s" << endl;
    if (duration.count() > max_time)
    {
        return size_results;
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


    // for (auto &it1 : kct)
    // {
    //     int u = it1.first;
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
    //         Weight u_v = cal_score(u_score, v_score, u_v_sim);
    //         // Weight u_v;
    //         // if (method == 0 || method == 1)
    //         // {
    //         //     u_v = cal_score(u_score, v_score, u_v_sim);
    //         // }
    //         // else if (method == 2)
    //         // {
    //         //     u_v = cal_score2(u_score, v_score, u_v_sim);
    //         // }

    //         if (!score_edges.count(u_v))
    //         {
    //             score_edges.emplace(u_v, vector<pair<int, int>>());
    //             // score_edges.emplace(u_v, unordered_map<int, set<int>>());
    //         }
    //         score_edges[u_v].push_back(pair<int, int>(u, v));
    //     }
    // }


    // cout << "calculate similarity complete. " << endl;
    // map<Weight, vector<vector<int>>, greater<Weight>> results;
    Weight ub;
    if (method == 0 || method == 1)
    {
        ub = Weight(2, 1);
    }
    else if (method == 2)
    {
        ub = Weight(2 * posKws.size(), 1);
    }
    // Weight ub(2, 1);
    map<Weight, vector<Edge>, greater<Weight>> pairsLeqUb;
    unordered_map<int, unordered_map<int, int>> S;
    vector<set<Edge>> nonDecSup; //[1] dont contain the old k-truss edges.
    nonDecSup.resize(2);
    auto it = s_truss_edges.begin();

    unordered_map<int, unordered_map<int, int>> T;
    vector<unordered_map<int, unordered_map<int, int>>> maxCTrusses;
    map<Edge, int> CTrussOfEdge;
    map<int, vector<int>> newPairs;
    queue<int> emptyIndexT;

    unordered_map<int, unordered_map<int, int>> updateEdgeGraph;

    unordered_map<int, unordered_map<int, int>> G;
    vector<vector<int>> maxCliques;                     // maximal cliques of G
    map<int, set<int>> cliquesOfNode;                   // the indices of maximal cliques containing each node in maxCliques
    set<int> notEmptyIndex, emptyIndexC, resultCliques; // the non-empty indices in maxCliques; the empty indices in maxCliques
    int cnt = 0;
    set<int> old_vertex;
    while (it != s_truss_edges.end() || !pairsLeqUb.empty())
    {
        // if (!G.empty())
        // {
        //     printGraphAsEdgeList(G);
        // }

        // for (auto &u_v : G)
        // {
        //     int u = u_v.first;
        //     for (auto &v_s : u_v.second)
        //     {
        //         int v = v_s.first;
        //         if (u == v)
        //         {
        //             int a = 1;
        //         }
        //     }
        // }

        // cout << "cnt: " << cnt << endl;
        // if (cnt == 82)
        // {
        //     int a = 1;
        // }
        cnt++;

        Weight w = it != s_truss_edges.end() ? it->first : WEIGHT_ZERO;

        vector<Edge> P_1;
        vector<Edge> P_2;
        vector<pair<Weight, int>> P_1_WeightIndex;
        if (it != s_truss_edges.end())
        {
            for (auto &e : it->second)
            {
                MaintainSupport_AddEdge(S, nonDecSup, e, k);
            }
            set<Edge> newEdges = GetMaxKTrussInc(T, k, nonDecSup);
            // if (w == Weight(5, 182))
            // {
            //     printGraphAsEdgeList(T);
            //     int a = 1;
            // }
            if (!newEdges.empty())
            {
                AddEdges(newEdges, T,
                         maxCTrusses,
                         CTrussOfEdge, newPairs,
                         emptyIndexT);

                set<int> nodeAttrs;
                // compute the pairwise weight for each pair in P
                for (auto &pair : newPairs)
                {
                    int u = pair.first;
                    Weight u_score = seq2att[u].second.first;
                    vector<int> u_att = seq2att[u].second.second;
                    auto &nodes = pair.second;
                    for (int v : nodes)
                    {
                        Weight v_score = seq2att[v].second.first;
                        vector<int> v_att = seq2att[v].second.second;
                        vector<int> com_att;
                        set_intersection(u_att.begin(), u_att.end(), v_att.begin(), v_att.end(), back_inserter(com_att));
                        Weight u_v_sim(com_att.size(), (u_att.size() + v_att.size() - com_att.size()));
                        Weight weight = cal_score(u_score, v_score, u_v_sim);
                        // Weight weight;
                        // if (method == 0 || method == 1)
                        // {
                        //     weight = cal_score(u_score, v_score, u_v_sim);
                        // }
                        // else if (method == 2)
                        // {
                        //     weight = cal_score2(u_score, v_score, u_v_sim);
                        // }
                        if (weight == WEIGHT_ZERO)
                        {
                            continue;
                        }

                        if (weight < ub)
                        {
                            pairsLeqUb[weight].push_back(make_pair(u, v));
                        }
                        else
                        {
                            P_2.push_back(make_pair(u, v));
                        }
                    }
                    nodeAttrs.clear();
                }
                newPairs.clear();
            }
            it++;
        }
        // cout << "add edges complete." << endl;
        // obtain P_1 and P_0(= P_0 \ P_1)
        for (auto weight_edges = pairsLeqUb.begin(); weight_edges != pairsLeqUb.end();)
        {
            Weight weight = weight_edges->first;
            if (weight >= w && weight < ub)
            {
                P_1.insert(P_1.end(), weight_edges->second.begin(), weight_edges->second.end());
                P_1_WeightIndex.push_back(make_pair(weight, weight_edges->second.size()));
                pairsLeqUb.erase(weight_edges++);
            }
            else
            {
                weight_edges++;
            }
        }

        for (auto &pair : P_2)
        {
            int u = pair.first;
            int v = pair.second;
            if (!updateEdgeGraph.count(u))
            {
                updateEdgeGraph.emplace(u, unordered_map<int, int>());
            }
            updateEdgeGraph[u][v] = 0;
            if (!updateEdgeGraph.count(v))
            {
                updateEdgeGraph.emplace(v, unordered_map<int, int>());
            }
            updateEdgeGraph[v][u] = 0;
        }
        // cout << "begin to niemch" << endl;
        AddEdges_NIEMCH(G, maxCliques, cliquesOfNode, notEmptyIndex, emptyIndexC,
                        resultCliques, updateEdgeGraph, start_time);
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - start_time;
        if (duration.count() > max_time)
        {
            return size_results;
        }
        // cout << "niemch complete" << endl;
        updateEdgeGraph.clear();
        vector<Edge>().swap(P_2);

        // lines 14-18 of Algorithm 6
        // QueryByPairs_IEMC(k, T, G, maxCliques, cliquesOfNode, notEmptyIndex, emptyIndexC,
        //                   resultCliques, P_1, P_1_WeightIndex, results, seq2att);

        int pairLIndex = 0;
        // unordered_map<int, unordered_map<int, int>> updateEdgeGraph;
        bool flag = false;
        for (auto &weight_cnt : P_1_WeightIndex)
        {
            Weight weight = weight_cnt.first;
            // cout << "weight: " << weight.toString() << endl;
            if (weight == Weight(5, 182))
            {
                int a = 1;
            }
            int cnt = weight_cnt.second;

            // 添加边到更新图中（边值初始化为-1）
            for (int i = pairLIndex; i < pairLIndex + cnt; i++)
            {
                int u = P_1[i].first;
                int v = P_1[i].second;
                updateEdgeGraph[u][v] = -1;
                updateEdgeGraph[v][u] = -1;
            }

            // 创建图的副本
            auto updateEdgeGraph1 = updateEdgeGraph;
            vector<vector<int>> newCliques;

            // newCliques = dynamicClique->AddEdges_NIEMCH(updateEdgeGraph);
            newCliques = AddEdges_NIEMCH(G, maxCliques, cliquesOfNode, notEmptyIndex, emptyIndexC,
                                         resultCliques, updateEdgeGraph, start_time);
            // cout << "update niemcn complete" << endl;
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = now - start_time;
            if (duration.count() > max_time)
            {
                return size_results;
            }
            // 获取结果
            map<Weight, vector<vector<int>>, greater<Weight>> results;
            GetResults(method, k, T, weight, newCliques, updateEdgeGraph1, results, seq2att);
            // cout << "get candidates complete" << endl;
            pairLIndex += cnt;

            // vector<set<int>> newVertex;
            // for (auto &w_r : results)
            // {

            //     for (auto &r : w_r.second)
            //     {
            //         set<int> new_vertex;
            //         sort(r.begin(), r.end());
            //         set_difference(r.begin(), r.end(), old_vertex.begin(), old_vertex.end(), inserter(new_vertex, new_vertex.begin()));
            //         old_vertex.insert(new_vertex.begin(), new_vertex.end());
            //         newVertex.push_back(new_vertex);
            //     }
            // }
            if (method == 0)
            {
                // map<Weight, vector<vector<int>>, greater<Weight>> size_results;
                for (auto &res : results)
                {
                    Weight score = res.first;
                    for (int i = 0; i < res.second.size(); i++)
                    {
                        // get subgraph
                        auto &clique = res.second[i];
                        vector<int> m = {3, 25, 168, 177};
                        bool f = true;
                        for (int i : m)
                        {
                            if (!count(clique.begin(), clique.end(), i))
                            {
                                f = false;
                                break;
                            }
                        }
                        if (f)
                        {
                            int a = 1;
                        }
                        unordered_map<int, unordered_map<int, int>> kct = GetSubGraph(T, clique);
                        kct = GetMaxKTruss2(kct, k);
                        unordered_set<int> visited;
                        DFS(kct, clique[0], visited);
                        if (kct.size() != clique.size() || visited.size() != clique.size())
                        {
                            cout << "error" << endl;
                        }
                        vector<int> M;
                        // for (auto &u_v : kct)
                        // {
                        //     int u = u_v.first;
                        //     for (auto &v_s : u_v.second)
                        //     {
                        //         int v = v_s.first;
                        //         if (u < v)
                        //             continue;
                        //         int sup = 0;
                        //         for (auto &w_s : kct[v])
                        //         {
                        //             int w = w_s.first;
                        //             if (u_v.second.count(w))
                        //             {
                        //                 sup++;
                        //             }
                        //         }
                        //         kct[u][v] = sup;
                        //         kct[v][u] = sup;
                        //     }
                        // }
                        // vector<int> m = {3, 6, 17, 18 ,19, 59, 71, 81, 177, 203 ,219, 304};
                        // bool f = true;
                        // for(int i : m){
                        //     if(!count(clique.begin(), clique.end(), i)){
                        //         f = false;
                        //         break;
                        //     }
                        // }
                        // if(f){
                        //     int a = 1;
                        // }
                        // if (count(clique.begin(), clique.end(), 177) && count(clique.begin(), clique.end(), 59) && count(clique.begin(), clique.end(), 1))
                        // {
                        //     int a = 1;
                        //     for(int u  : clique){
                        //         cout << u << ", ";
                        //     }
                        //     cout << endl;
                        //     printGraphAsEdgeList(kct);
                        // }
                        // NaiveEnum(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        // NaiveEnum_Truss(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        // NaiveEnum_CntTruss(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        NaiveEnum_CntTrussNsize(score, M, clique, kct, seq2att, k, sthd1, sthd2, size_results, start_time);
                        // NaiveEnum_CntTrussNsize(newVertex[i], score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                        // NaiveEnum_CntTrussNsizeVert(score, M, clique, kct, seq2att, k, opt, sthd1, sthd2, size_results);
                    }
                }
                if (!size_results.empty())
                {
                    flag = true;
                    printf_results(dg, size_results, seq2att);
                    break;
                }
            }
            else
            {
                if (!results.empty())
                {

                    flag = true;
                    printf_results(dg, results, seq2att);
                    size_results = results;
                    break;
                }
            }
            results.clear();

            // 清空图
            updateEdgeGraph.clear();
            updateEdgeGraph1.clear();
        }

        ub = w;

        // resultCnt += results.size();
        // if (resultCnt > 10)
        // {
        //     break;
        // }
        if (flag)
        {
            // printf_results(size_results, seq2att);
            break;
        }
        // results.clear();
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // cout << "search time: " << duration.count() << "s" << endl;
    // for (auto &res : results)
    // {
    //     cout << "Weight: " << res.first.toString() << ", ";
    //     for (auto &clique : res.second)
    //     {
    //         for (auto &node : clique)
    //         {
    //             cout << node << " ";
    //         }
    //         cout << endl;
    //     }
    // }
    return move(size_results);
}