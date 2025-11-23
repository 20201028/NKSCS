#include "KScore.h"
void reorderEL(vector<Edge> &sorted_elbys, map<Edge, int> &sorted_ep, map<Edge, Weight> &supd, unordered_map<int, int> &svp, Edge &e1, map<Weight, int> &w2id, Weight &cur_s)
{

    Weight val = supd[e1];
    int pos1 = sorted_ep[e1];
    int id = w2id[val];
    int cp = svp[id];
    if (cp != pos1)
    {
        // 把开始的边与该边交换
        Edge tmp2 = sorted_elbys[cp];
        sorted_ep[e1] = cp;
        sorted_elbys[cp] = e1;
        sorted_ep[tmp2] = pos1;
        sorted_elbys[pos1] = tmp2;
        svp[id] = cp + 1;
    }
    else
    {
        if (sorted_elbys.size() > cp + 1 && supd[sorted_elbys[cp + 1]] == val) // 是第一条边
        {
            svp[id] = cp + 1;
        }
        else // 该值下仅仅这条边
        {
            svp[id] = -1;
        }
    }
    // int cur_s_id = w2id[cur_s];
    // if (svp.find(cur_s_id) == svp.end() || svp[cur_s_id] == -1)
    //     svp[cur_s_id] = cp;
    supd[e1] = cur_s;
}
// KScoreTrussDecomposition::KScoreTrussDecomposition(unordered_map<int, unordered_map<int, int>> &graph, map<Edge, Weight> graph_sim, int k) // edge's first vertex less than second
// {
//     // Weight s_max = WEIGHT_ZERO;

//     // for (auto &it : graph_sim)
//     // {
//     //     Weight s = it.second;
//     //     if(!w2id.count(s)){
//     //         id2w.emplace(id, s);
//     //         w2id.emplace(s, id);
//     //         id++;
//     //     }
//     //     if (s > s_max)
//     //     {
//     //         s_max = s;
//     //     }
//     // }
//     set<Weight> scores;
//     vector<pair<Edge, int>> changed_edges;
//     for (auto &it : graph_sim)
//     {
//         Weight s = it.second;
//         scores.insert(s);
//     }
//     // sort(scores.begin(), scores.end(), [](const Weight &a, const Weight &b)
//     //      { return a < b; });
//     map<int, Weight> id2w;
//     map<Weight, int> w2id;
//     int id = 0;
//     for (auto &it : scores)
//     {

//         id2w.emplace(id, it);
//         w2id.emplace(it, id);
//         id++;
//     }
//     // bucket Sort edge
//     vector<Edge> sorted_elbys(graph_sim.size()); // sorted_elbys和sorted_ep相互对应
//     map<Edge, int> sorted_ep;                    // 边在sorted_elbys中的位置
//     unordered_map<int, int> svp;
//     vector<int> bucket;
//     int s_num = w2id.size();
//     bucket.resize(s_num, 0);
//     for (auto &e : graph_sim)
//     {
//         int id = w2id[e.second];
//         bucket[id]++; // bucket[i] = the number of edges whose support = i
//     }
//     int p = 0;
//     for (int j = 0; j < s_num; j++)
//     {
//         int tmp = bucket[j];
//         bucket[j] = p; // record the first index of edge whose support=j after sort
//         p += tmp;
//     }
//     for (auto &e : graph_sim)
//     {
//         int id = w2id[e.second];
//         sorted_elbys[bucket[id]] = e.first;     // edges sort as support increasing
//         sorted_ep.emplace(e.first, bucket[id]); // edge : index after sort
//         svp.emplace(id, bucket[id]);            // support: the first index of edge whose support=support after sort
//         bucket[id]++;
//     }
//     //
//     // Weight cur_s = *scores.begin();
//     // set<Edge> kedgelist;
//     // map<int, set<Edge>> k2edge;
//     for (int i = 0; i < sorted_elbys.size(); i++)
//     {
//         Edge e = sorted_elbys[i];
//         Weight cur_s = graph_sim[e];
//         // if (val > cur_s)
//         // {
//         //     // k2edge.emplace(kMax, kedgelist);
//         //     cur_s = val + 2;
//         //     // kedgelist.clear();
//         // }
//         int src = e.first;
//         int dst = e.second;
//         if (graph[src].size() > graph[dst].size())
//         {
//             swap(src, dst);
//         }
//         for (auto it : graph[src])
//         {
//             int w = it.first;
//             if (graph[dst].count(w))
//             {
//                 // update support, if a edge be less than k, then change its sim be cur_s
//                 Edge e1;
//                 if (src < w)
//                 {
//                     e1 = make_pair(src, w);
//                 }
//                 else
//                 {
//                     e1 = make_pair(w, src);
//                 }
//                 if (graph_sim[e1] > cur_s && graph[src][w] <= k - 2)
//                 {
//                     reorderEL(sorted_elbys, sorted_ep, graph_sim, svp, e1, w2id, cur_s);
//                 }
//                 Edge e2;
//                 if (dst < w)
//                 {
//                     e2 = make_pair(dst, w);
//                 }
//                 else
//                 {
//                     e2 = make_pair(w, dst);
//                 }
//                 if (graph_sim[e2] > cur_s && graph[dst][w] <= k - 2)
//                 {

//                     reorderEL(sorted_elbys, sorted_ep, graph_sim, svp, e2, w2id, cur_s);
//                 }
//                 graph[src][w]--;
//                 graph[w][src]--;
//                 graph[dst][w]--;
//                 graph[w][dst]--;
//             }
//         }
//         changed_edges.push_back(make_pair(e, graph[src][dst]));
//         graph[src].erase(dst);
//         graph[dst].erase(src);

//         // kedgelist.insert(e);
//         edge2score.emplace(e, cur_s);
//     }
//     while (!changed_edges.empty())
//     {
//         auto &it = changed_edges.back();
//         Edge e = it.first;
//         int sup = it.second;
//         graph[e.first][e.second] = sup;
//         graph[e.second][e.first] = sup;
//         changed_edges.pop_back();
//     }
//     // k2edge.emplace(kMax, kedgelist);
// }
void update_sim(vector<Edge> &edges, map<Weight, vector<Edge>> &sim_graph, map<Edge, Weight> &graph_sim, Edge &e, Weight &cur_s)
{
    Weight s = graph_sim[e];
    sim_graph[s].erase(find(sim_graph[s].begin(), sim_graph[s].end(), e));
    graph_sim[e] = cur_s;
    edges.push_back(e);
}
KScoreTrussDecomposition::KScoreTrussDecomposition(unordered_map<int, unordered_map<int, int>> &graph, map<Weight, vector<Edge>> &sim_graph, map<Edge, Weight> &graph_sim, int k,
                                                   chrono::high_resolution_clock::time_point startTime)
{
    // vector<pair<Edge, int>> changed_edges;
    for (auto &it : sim_graph)
    {
        Weight cur_s = it.first;
        vector<Edge> &edges = it.second;
        while (!edges.empty())
        {
            Edge e = edges.back();

            int src = e.first;
            int dst = e.second;
            if (src == 17 && dst == 203)
            {
                int a = 1;
            }
            edges.pop_back();
            if (graph[src].size() > graph[dst].size())
            {
                swap(src, dst);
            }
            for (auto it : graph[src])
            {
                auto now = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = now - startTime;
                if (duration.count() > max_time)
                {
                    edge2score.clear();
                    return;
                }
                int w = it.first;
                if (graph[dst].count(w))
                {
                    // update support, if a edge be less than k, then change its sim be cur_s
                    Edge e1;
                    if (src < w)
                    {
                        e1 = make_pair(src, w);
                    }
                    else
                    {
                        e1 = make_pair(w, src);
                    }
                    if (graph_sim[e1] > cur_s && graph[src][w] == k - 2)
                    {
                        update_sim(edges, sim_graph, graph_sim, e1, cur_s);
                    }
                    Edge e2;
                    if (dst < w)
                    {
                        e2 = make_pair(dst, w);
                    }
                    else
                    {
                        e2 = make_pair(w, dst);
                    }
                    if (graph_sim[e2] > cur_s && graph[dst][w] == k - 2)
                    {

                        update_sim(edges, sim_graph, graph_sim, e2, cur_s);
                    }
                    graph[src][w]--;
                    graph[w][src]--;
                    graph[dst][w]--;
                    graph[w][dst]--;
                }
            }
            // changed_edges.push_back(make_pair(e, graph[src][dst]));
            graph[src].erase(dst);
            graph[dst].erase(src);

            // kedgelist.insert(e);
            edge2score.emplace(e, cur_s);
        }
    }
    // while (!changed_edges.empty())
    // {
    //     auto &it = changed_edges.back();
    //     Edge e = it.first;
    //     int sup = it.second;
    //     graph[e.first][e.second] = sup;
    //     graph[e.second][e.first] = sup;
    //     changed_edges.pop_back();
    // }
}