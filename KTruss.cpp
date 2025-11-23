#include "KTruss.h"
void reorderEL(vector<Edge> &sorted_elbys, map<Edge, int> &sorted_ep, map<Edge, int> &supd, unordered_map<int, int> &svp, Edge &e1)
{
    int val = supd[e1];
    int pos1 = sorted_ep[e1];
    int cp = svp[val];
    if (cp != pos1)
    {
        Edge tmp2 = sorted_elbys[cp];
        sorted_ep[e1] = cp;
        sorted_elbys[cp] = e1;
        sorted_ep[tmp2] = pos1;
        sorted_elbys[pos1] = tmp2;
        svp[val] = cp + 1;
    }
    else
    {
        if (sorted_elbys.size() > cp + 1 && supd[sorted_elbys[cp + 1]] == val)
        {
            svp[val] = cp + 1;
        }
        else
        {
            svp[val] = -1;
        }
    }
    if (svp.find(val - 1) == svp.end() || svp[val - 1] == -1)
        svp[val - 1] = cp;
    supd[e1] = val - 1;
}
TrussDecomposition::TrussDecomposition() { kMax = 2; }
TrussDecomposition::TrussDecomposition(DataGraph datagraph)
{
    if (datagraph.id2seq.empty())
    {
        kMax = 0;
        return;
    }
    // total_order.clear();
    // degree.clear();
    map<Edge, int> support;
    set<Edge> visited;
    // trussd.clear();
    // for (int i = 0; i < datagraph.AdjList.size(); i++)
    // {
    //     total_order.push_back(i);
    //     int deg = datagraph.AdjList[i].size();
    //     degree.push_back(deg);
    // }
    // sort(total_order.begin(), total_order.end(), [this](int a, int b)
    //      { return degree[a] > degree[b] || (degree[a] == degree[b] && a > b); });
    // calculate support
    // int n = total_order.size(), nBins = 0, mp = 0;
    int kmax = 0;
    for (int v = 0; v < datagraph.AdjList.size(); v++)
    {
        vector<int> vNei = datagraph.AdjList[v];
        sort(vNei.begin(), vNei.end());
        for (int j = 0; j < datagraph.AdjList[v].size(); j++)
        {
            int u = datagraph.AdjList[v][j];
            if (v > u)
                continue;
            vector<int> uNei = datagraph.AdjList[u];
            sort(uNei.begin(), uNei.end());
            vector<int> com;
            set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), back_inserter(com));
            int sup = com.size();
            Edge e = make_pair(v, u);
            support.emplace(e, sup);
            if (sup > kmax)
                kmax = sup;
        }
    }
    // sup = support;
    // bucket Sort edge
    vector<Edge> sorted_elbys(support.size());
    map<Edge, int> sorted_ep;
    unordered_map<int, int> svp;
    vector<int> bucket;

    bucket.resize(kmax + 1, 0);
    for (auto &e : support)
    {
        bucket[e.second]++; // bucket[i] = the number of edges whose support = i, 同一支持度的边有多少
    }
    int p = 0;
    for (int j = 0; j < kmax + 1; j++)
    {
        int tmp = bucket[j];
        bucket[j] = p; // record the first index of edge whose support=j after sort
        p += tmp;
    }
    for (auto &e : support)
    {
        sorted_elbys[bucket[e.second]] = e.first;     // edges sort as support increasing
        sorted_ep.emplace(e.first, bucket[e.second]); // edge : index after sort
        svp.emplace(e.second, bucket[e.second]);      // support: the first index of edge whose support=support after sort
        bucket[e.second]++;
    }
    //
    kMax = 2;
    vector<Edge> kedgelist;
    // map<int, set<Edge>> k2edge;
    for (int i = 0; i < sorted_elbys.size(); i++)
    {
        Edge e = sorted_elbys[i];
        int val = support[e];
        if (val > kMax - 2)
        {
            k2edge.emplace(kMax, kedgelist);
            kMax = val + 2;
            kedgelist.clear();
        }
        int src = e.first;
        int dst = e.second;
        vector<int> srcNei = datagraph.AdjList[src];
        sort(srcNei.begin(), srcNei.end());
        vector<int> dstNei = datagraph.AdjList[dst];
        sort(dstNei.begin(), dstNei.end());
        vector<int> com;
        set_intersection(srcNei.begin(), srcNei.end(), dstNei.begin(), dstNei.end(), back_inserter(com));
        for (auto &w : com)
        {
            Edge e1;
            if (w < src)
                e1 = make_pair(w, src);
            else
                e1 = make_pair(src, w);

            Edge e2;
            if (w < dst)
                e2 = make_pair(w, dst);
            else
                e2 = make_pair(dst, w);
            if (visited.find(e1) == visited.end() && visited.find(e2) == visited.end())
            {
                if (support[e1] > kMax - 2)
                    reorderEL(sorted_elbys, sorted_ep, support, svp, e1);
                if (support[e2] > kMax - 2)
                    reorderEL(sorted_elbys, sorted_ep, support, svp, e2);
            }
        }
        kedgelist.push_back(e);
        visited.insert(e);
    }
    k2edge.emplace(kMax, kedgelist);
}
// TrussDecomposition::TrussDecomposition(unordered_map<int, unordered_map<int, int>> &graph, int k)
// {
//     map<Edge, int> support;
//     int kmax = 0;
//     for (auto it : graph)
//     {
//         int u = it.first;
//         for (auto it2 : it.second)
//         {
//             int v = it2.first;
//             if (u > v)
//                 continue;
//             Edge e = make_pair(u, v);
//             int sup = it2.second;
//             support.emplace(e, sup);
//             if (sup > kmax)
//                 kmax = sup;
//         }
//     }
//     // for (auto &it1 : graph)
//     // {
//     //     int v = it1.first;
//     //     for (auto &it2 : it1.second)
//     //     {
//     //         int u = it2.first;
//     //         if (v > u)
//     //             continue;
//     //         int sup = 0;
//     //         for(auto &it3 : graph[u]){
//     //             int w = it3.first;
//     //             if(graph[v].count(w)){
//     //                 sup++;
//     //             }
//     //         }
//     //         Edge e = make_pair(v, u);
//     //         support.emplace(e, sup);
//     //         if (sup > kmax)
//     //             kmax = sup;
//     //     }
//     // }
//     // sup = support;
//     // bucket Sort edge
//     vector<Edge> sorted_elbys(support.size());
//     map<Edge, int> sorted_ep;
//     unordered_map<int, int> svp;
//     vector<int> bucket;

//     bucket.resize(kmax + 1, 0);
//     for (auto &e : support)
//     {
//         bucket[e.second]++; // bucket[i] = the number of edges whose support = i
//     }
//     int p = 0;
//     for (int j = 0; j < kmax + 1; j++)
//     {
//         int tmp = bucket[j];
//         bucket[j] = p; // record the first index of edge whose support=j after sort
//         p += tmp;
//     }
//     for (auto &e : support)
//     {
//         sorted_elbys[bucket[e.second]] = e.first;     // edges sort as support increasing
//         sorted_ep.emplace(e.first, bucket[e.second]); // edge : index after sort
//         svp.emplace(e.second, bucket[e.second]);      // support: the first index of edge whose support=support after sort
//         bucket[e.second]++;
//     }
//     //
//     kMax = 2;
//     // set<Edge> kedgelist;
//     // map<int, set<Edge>> k2edge;
//     for (int i = 0; i < sorted_elbys.size(); i++)
//     {
//         Edge e = sorted_elbys[i];
//         int val = support[e];
//         if (val > kMax - 2)
//         {
//             // k2edge.emplace(kMax, kedgelist);
//             kMax = val + 2;
//             // kedgelist.clear();
//         }
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
//                 Edge e1;
//                 if (w < src)
//                     e1 = make_pair(w, src);
//                 else
//                     e1 = make_pair(src, w);

//                 Edge e2;
//                 if (w < dst)
//                     e2 = make_pair(w, dst);
//                 else
//                     e2 = make_pair(dst, w);
//                 if (!(edge_truss.count(w) && edge_truss[w].count(dst)) && !(edge_truss.count(w) && edge_truss[w].count(src)))
//                 {
//                     if (support[e1] > kMax - 2)
//                         reorderEL(sorted_elbys, sorted_ep, support, svp, e1);
//                     if (support[e2] > kMax - 2)
//                         reorderEL(sorted_elbys, sorted_ep, support, svp, e2);
//                 }
//             }
//         }

//         // kedgelist.insert(e);
//         // trussd.emplace(e, kMax);
//         auto it = edge_truss.find(src);
//         if(it == edge_truss.end()){
//             edge_truss.emplace(src, unordered_map<int, int>(dst, kMax));
//         }
//         else{
//             it->second.emplace(dst, kMax);
//         }
//         auto it2 = edge_truss.find(dst);
//         if(it2 == edge_truss.end()){
//             edge_truss.emplace(dst, unordered_map<int, int>(src, kMax));
//         }
//         else{
//             it2->second.emplace(src, kMax);
//         }
//     }
//     // k2edge.emplace(kMax, kedgelist);
// }
TrussDecomposition::TrussDecomposition(DataGraph datagraph, int k, chrono::high_resolution_clock::time_point startTime)
{
    // total_order.clear();
    // degree.clear();
    map<Edge, int> support;
    int sup_max = 0;
    for (int v = 0; v < datagraph.AdjList.size(); v++)
    {
        vector<int> vNei = datagraph.AdjList[v];
        sort(vNei.begin(), vNei.end());
        for (int j = 0; j < datagraph.AdjList[v].size(); j++)
        {
            int u = datagraph.AdjList[v][j];
            if (v > u)
                continue;
            vector<int> uNei = datagraph.AdjList[u];
            sort(uNei.begin(), uNei.end());
            vector<int> com;
            set_intersection(vNei.begin(), vNei.end(), uNei.begin(), uNei.end(), back_inserter(com));
            int sup = com.size();
            Edge e = make_pair(v, u);
            support.emplace(e, sup);
            if (sup > sup_max)
                sup_max = sup;
        }
    }
    // sup = support;
    // bucket Sort edge
    vector<Edge> sorted_elbys(support.size());
    map<Edge, int> sorted_ep;
    unordered_map<int, int> svp;
    vector<int> bucket;

    bucket.resize(sup_max + 1, 0);
    for (auto &e : support)
    {
        bucket[e.second]++; // bucket[i] = the number of edges whose support = i
    }
    int p = 0;
    for (int j = 0; j < sup_max + 1; j++)
    {
        int tmp = bucket[j];
        bucket[j] = p; // record the first index of edge whose support=j after sort
        p += tmp;
    }
    for (auto &e : support)
    {
        sorted_elbys[bucket[e.second]] = e.first;     // edges sort as support increasing
        sorted_ep.emplace(e.first, bucket[e.second]); // edge : index after sort
        svp.emplace(e.second, bucket[e.second]);      // support: the first index of edge whose support=support after sort
        bucket[e.second]++;
    }
    vector<Edge> visited;
    // map<int, set<Edge>> k2edge;
    int kIndex = 0;
    for (int i = 0; i < sorted_elbys.size(); i++)
    {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = now - startTime;
        if (duration.count() > max_time)
        {
            cout << "Truss Time out" << endl;
            kTruss.clear();
            return;
        }
        Edge e = sorted_elbys[i];
        int val = support[e];
        if (val < k - 2)
        {
            int src = e.first;
            int dst = e.second;
            vector<int> srcNei = datagraph.AdjList[src];
            sort(srcNei.begin(), srcNei.end());
            vector<int> dstNei = datagraph.AdjList[dst];
            sort(dstNei.begin(), dstNei.end());
            vector<int> com;
            set_intersection(srcNei.begin(), srcNei.end(), dstNei.begin(), dstNei.end(), back_inserter(com));
            for (auto &w : com)
            {
                Edge e1;
                if (w < src)
                    e1 = make_pair(w, src);
                else
                    e1 = make_pair(src, w);

                Edge e2;
                if (w < dst)
                    e2 = make_pair(w, dst);
                else
                    e2 = make_pair(dst, w);
                if (find(visited.begin(), visited.end(), e1) == visited.end() && find(visited.begin(), visited.end(), e2) == visited.end())
                {
                    if (support[e1] >= k - 2)
                        reorderEL(sorted_elbys, sorted_ep, support, svp, e1);
                    if (support[e2] >= k - 2)
                        reorderEL(sorted_elbys, sorted_ep, support, svp, e2);
                }
            }
            visited.push_back(e);
        }
        else
        {
            // int u_id = datagraph.seq2id[e.first];
            // int v_id = datagraph.seq2id[e.second];
            // if (!kTruss.count(u_id))
            // {
            //     kTruss.emplace(u_id, unordered_map<int, int>());
            // }
            // kTruss[u_id].emplace(v_id, support[e]);
            // if (kTruss.count(v_id) == 0)
            // {
            //     kTruss.emplace(v_id, unordered_map<int, int>());
            // }
            // kTruss[v_id].emplace(u_id, support[e]);
            if (!kTruss.count(e.first))
            {
                kTruss.emplace(e.first, unordered_map<int, int>());
            }
            kTruss[e.first].emplace(e.second, support[e]);
            if (kTruss.count(e.second))
            {
                kTruss.emplace(e.second, unordered_map<int, int>());
            }
            kTruss[e.second].emplace(e.first, support[e]);
        }
    }
}

// void TrussDecomposition::UpdateTrussofInsert2(vector<Edge> &newEs, Edge &e, DataGraph *relSubG)
// {
//     // if (e.first == 1856 && e.second == 33749)
//     // {
//     //     cout << "here" << endl;
//     // }
//     // k1,k2
//     int v = relSubG->id2seq[e.first];
//     int u = relSubG->id2seq[e.second];

//     if (relSubG->AdjList[u].size() > relSubG->AdjList[v].size())
//         swap(u, v);
//     vector<int> vN = relSubG->AdjList[v];
//     vector<int> uN = relSubG->AdjList[u];
//     vector<int> K(this->kMax + 2);
//     int maxTruss = this->kMax + 1;
//     for (int w : uN)
//     {
//         if (count(vN.begin(), vN.end(), w))
//         {
//             // (v,w) and (u,w)
//             int t1 = this->trussd[make_pair(min(v, w), max(v, w))];
//             maxTruss = max(maxTruss, t1);
//             int t2 = this->trussd[make_pair(min(u, w), max(u, w))];
//             maxTruss = max(maxTruss, t2);
//             int mint = min(t1, t2);
//             for (int i = 2; i <= mint; i++)
//                 K[i]++;
//         }
//     }
//     int K1 = 2, K2 = 2;
//     for (int i = 2; i <= maxTruss; i++)
//     {
//         if (K[i] >= i - 2)
//             K1 = i;
//         if (K[i - 1] >= i - 2)
//             K2 = i;
//     }
//     //////////////////////////////////////////////
//     int K_max = K2 - 1;

//     //  cout << "kmax= " << K_max <<endl;

//     map<int, vector<Edge>> L;
//     Edge e_uv = make_pair(min(u, v), max(u, v));
//     this->trussd.emplace(e_uv, K1);
//     if (K1 > this->kMax)
//         this->kMax = K1;
//     newEs.push_back(e_uv);
//     // if (K1 == K_max)
//     // {
//     //     vector<Edge> t = {e_uv};
//     //     L.emplace(K_max, t);
//     // }
//     // if (relSubG->AdjList[u].size() > relSubG->AdjList[v].size())
//     //     swap(u, v);
//     // vN = relSubG->AdjList[v];
//     // uN = relSubG->AdjList[u];
//     for (int w : uN)
//     {
//         if (count(vN.begin(), vN.end(), w))
//         {
//             // (v,w) and (u,w)
//             Edge e1 = make_pair(min(v, w), max(v, w));
//             int t1 = this->trussd[e1];
//             Edge e2 = make_pair(min(u, w), max(u, w));
//             int t2 = this->trussd[e2];
//             int mint = min(t1, t2);
//             if (mint <= K_max)
//             {
//                 auto it = L.find(mint);
//                 if (t1 == mint)
//                 {

//                     if (it == L.end())
//                     {
//                         vector<Edge> t = {e1};
//                         L.emplace(mint, t);
//                     }
//                     else
//                     {
//                         it->second.push_back(e1);
//                     }
//                 }

//                 if (t2 == mint)
//                 {
//                     if (it == L.end())
//                     {
//                         vector<Edge> t = {e2};
//                         L.emplace(mint, t);
//                     }
//                     else
//                     {
//                         it->second.push_back(e2);
//                     }
//                 }
//             }
//         }
//     }

//     //  cout << "kmax= " << K_max <<endl;
//     map<Edge, int> _sup;
//     for (int k = K_max; k >= 2; k--)
//     {
//         //     cout << "L[k].size= " << L[k].size() <<endl;

//         vector<Edge> &kE = L[k];
//         queue<Edge> q1;
//         for (const auto &edge : kE)
//         {
//             q1.push(edge);
//         }
//         while (!q1.empty())
//         {
//             Edge e_xy = q1.front();
//             q1.pop();
//             _sup.emplace(e_xy, 0);
//             int u = e_xy.first, v = e_xy.second;
//             if (relSubG->AdjList[u].size() > relSubG->AdjList[v].size())
//                 swap(u, v);
//             vN = relSubG->AdjList[v];
//             uN = relSubG->AdjList[u];
//             for (int w : uN)
//             {
//                 if (count(vN.begin(), vN.end(), w))
//                 {
//                     // (v,w) and (u,w)
//                     Edge e1 = make_pair(min(v, w), max(v, w));
//                     int t1 = this->trussd[e1];
//                     Edge e2 = make_pair(min(u, w), max(u, w));
//                     int t2 = this->trussd[e2];
//                     if (t1 < k || t2 < k)
//                         continue;
//                     _sup[e_xy]++;

//                     if (t1 == k && find(kE.begin(), kE.end(), e1) == kE.end())
//                     {
//                         q1.push(e1);
//                         kE.push_back(e1);
//                     }

//                     if (t2 == k && find(kE.begin(), kE.end(), e2) == kE.end())
//                     {
//                         q1.push(e2);
//                         kE.push_back(e2);
//                     }
//                 }
//             }
//         }
//         // vector<Edge> q = kE;
//         for (const auto &edge : kE)
//         {
//             q1.push(edge);
//         }
//         // auto qit = q.begin();
//         // while (qit != q.end())
//         while (!q1.empty())
//         {

//             // Edge e_xy = *qit;
//             Edge e_xy = q1.front();
//             q1.pop();
//             if (_sup[e_xy] <= k - 2)
//             {
//                 int u = e_xy.first, v = e_xy.second;
//                 // if (e.first == 34 && e.second == 57)
//                 // {
//                 //     cout << "here" << endl;
//                 // }
//                 auto it = find(kE.begin(), kE.end(), e_xy);
//                 if (it != kE.end())
//                 {
//                     // 删除找到的元素
//                     kE.erase(it);
//                 }
//                 if (relSubG->AdjList[u].size() > relSubG->AdjList[v].size())
//                     swap(u, v);
//                 vN = relSubG->AdjList[v];
//                 uN = relSubG->AdjList[u];
//                 for (int w : uN)
//                 {
//                     if (count(vN.begin(), vN.end(), w))
//                     {
//                         // (v,w) and (u,w)
//                         Edge e1 = make_pair(min(v, w), max(v, w));
//                         int t1 = this->trussd[e1];
//                         Edge e2 = make_pair(min(u, w), max(u, w));
//                         int t2 = this->trussd[e2];
//                         if (t1 < k || t2 < k)
//                             continue;
//                         auto it1 = find(kE.begin(), kE.end(), e1);
//                         auto it2 = find(kE.begin(), kE.end(), e2);
//                         if (t1 == k && it1 == kE.end())
//                             continue;
//                         if (t2 == k && it2 == kE.end())
//                             continue;
//                         if (it1 != kE.end())
//                         {
//                             _sup[e1]--;
//                             if (_sup[e1] <= k - 2)
//                             {
//                                 // auto qit1 = find(q.begin(), q.end(), e1);
//                                 // if (qit1 < qit)
//                                 //     q.push_back(e1);
//                                 q1.push(e1);
//                             }
//                         }
//                         if (it2 != kE.end())
//                         {
//                             _sup[e2]--;
//                             if (_sup[e2] <= k - 2)
//                             {
//                                 // auto qit2 = find(q.begin(), q.end(), e2);
//                                 // if (qit2 < qit)
//                                 //     q.push_back(e2);
//                                 q1.push(e2);
//                             }
//                         }
//                     }
//                 }
//             }

//             // qit++;
//         }
//         for (auto &e_xy : kE)
//         {
//             int kt = ++this->trussd[e_xy];
//             if (kt > this->kMax)
//                 this->kMax = kt;
//         }
//     }

//     // Finally, update the truss value for all remaining edges in L_k
// }