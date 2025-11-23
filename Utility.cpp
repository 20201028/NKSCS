#include "Utility.h"
#include <stack>
/*method 0: PN + NN, 1: PN */
unordered_map<int, pair<bool, pair<Weight, vector<int>>>> LoadAtt(int method, string AttFileName, vector<int> posKws, vector<int> negKws)
{
    unordered_map<int, pair<bool, pair<Weight, vector<int>>>> id2att;
    // string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    ifstream ain(AttFileName.c_str());

    if (!ain)
    {
        cout << "Fail to read " << AttFileName << "." << endl;
        return {};
    }

    string aline;

    while (getline(ain, aline))
    {
        if (aline.find('#') != string::npos)
            continue;
        string ver_s = aline.substr(0, aline.find("\t"));
        string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
        int ver = stoi(ver_s); // string2float
        int start = 0, end = 0;
        vector<int> att;
        int w = 0;
        while ((end = att_s.find(",", start)) != string::npos)
        {
            w = stoi(att_s.substr(start, end - start));
            att.push_back(w);
            start = end + 1;
        }
        w = stoi(att_s.substr(start));
        att.push_back(w);
        vector<int> comKw;
        set_intersection(posKws.begin(), posKws.end(), att.begin(), att.end(), back_inserter(comKw));
        if (comKw.size() > 0)
        {
            int PN = comKw.size();
            comKw.clear();
            set_intersection(negKws.begin(), negKws.end(), att.begin(), att.end(), back_inserter(comKw));
            // pair<Weight, vector<int>> kscoreAtt(Weight(PN, NN + posKws.size()), att);
            pair<Weight, vector<int>> kscoreAtt;
            int NN = comKw.size();
            if (method == 0)
            {
                kscoreAtt.first = Weight(PN, NN + posKws.size());
                kscoreAtt.second = att;
            }
            else
            {
                kscoreAtt.first = Weight(PN, posKws.size());
                kscoreAtt.second = att;
            }
            // pair<Weight, vector<int>> kscoreAtt(Weight(PN, NN + 1), att);
            // pair<double, vector<int>> kscoreAtt(PN, att);
            if (NN == 0)
            {
                // id2att.emplace(ver, pair<bool, pair<Weight, vector<int>>>(false, kscoreAtt));
                id2att.emplace(ver, pair<bool, pair<Weight, vector<int>>>(true, kscoreAtt));
            }
            else
            {
                // id2att.emplace(ver, pair<bool, pair<Weight, vector<int>>>(true, kscoreAtt));
                id2att.emplace(ver, pair<bool, pair<Weight, vector<int>>>(false, kscoreAtt));
            }
        }
    }
    ain.close();
    cout << "Loaded attribute successfully!" << endl;
    return move(id2att);
}
DataGraph LoadGraph(string dataset, unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &id2att)
{
    DataGraph datagraph;
    string StructFileName = "DataGraph/" + dataset + "/" + "graph.txt";
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
        string src_s = sline.substr(0, sline.find("\t"));
        string dst_s = sline.substr(sline.find("\t") + 1, sline.find("\n") - sline.find("\t") - 1);
        int src = stoi(src_s); // string2float
        int dst = stoi(dst_s);
        if (id2att.count(src) && id2att.count(dst))
        {
            datagraph.addEdgeNoMatinC(src, dst);
        }
    }

    cout << "Loaded dataset successfully!" << endl;
    return move(datagraph);
}
Weight cal_score(Weight u_score, Weight v_score, Weight u_v_sim)
{
    // return move(u_score * v_score * u_v_sim);
    return move((u_score + v_score) * u_v_sim);
}
void DFS(const unordered_map<int, unordered_map<int, int>> &kct, int start, unordered_set<int> &visited, vector<int> &component)
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
            component.push_back(node);
            for (const auto &neighbor : kct.at(node))
            {
                stack.push(neighbor.first);
            }
        }
    }
}
vector<vector<int>> FindConnectedComponents(const unordered_map<int, unordered_map<int, int>> &kct)
{
    vector<vector<int>> components;
    unordered_set<int> visited;

    for (const auto &node : kct)
    {
        if (visited.find(node.first) == visited.end())
        {
            vector<int> component;
            DFS(kct, node.first, visited, component);
            components.push_back(component);
        }
    }

    return move(components);
}
void printGraphAsEdgeList(const std::unordered_map<int, std::unordered_map<int, int>> &graph)
{
    std::cout << "edges = [\n";
    for (auto &u_v : graph)
    {
        int src = u_v.first;
        for (const auto &[dst, weight] : graph.at(src))
        {
            std::cout << "    (" << src << ", " << dst << ", " << weight << "),\n";
        }
    }
    std::cout << "]\n";
}
unordered_map<int, vector<int>> LoadAtt(string dataset)
{
    unordered_map<int, vector<int>> id2att;
    string AttFileName = "DataGraph/" + dataset + "/" + "user_attributes_int.txt";
    ifstream ain(AttFileName.c_str());

    if (!ain)
    {
        cout << "Fail to read " << AttFileName << "." << endl;
        return {};
    }

    string aline;

    while (getline(ain, aline))
    {
        if (aline.find('#') != string::npos)
            continue;
        string ver_s = aline.substr(0, aline.find("\t"));
        string att_s = aline.substr(aline.find("\t") + 1, aline.find("\n") - aline.find("\t") - 1);
        int ver = stoi(ver_s); // string2float
        int start = 0, end = 0;
        vector<int> att;
        int w = 0;
        while ((end = att_s.find(",", start)) != string::npos)
        {
            w = stoi(att_s.substr(start, end - start));
            att.push_back(w);
            start = end + 1;
        }
        w = stoi(att_s.substr(start));
        att.push_back(w);
        id2att.emplace(ver, att);
    }
    ain.close();
    cout << "Loaded attribute successfully!" << endl;
    return move(id2att);
}
void calculate_score_edges(unordered_map<int, unordered_map<int, int>> &kct,
                           unordered_map<int, pair<bool, pair<Weight, vector<int>>>> &seq2att,
                           map<Edge, Weight> &edge_sim,
                           map<Weight, vector<pair<int, int>>> &score_edges)
{
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
            Edge e(u, v);
            score_edges[u_v].push_back(e);
            edge_sim.emplace(e, u_v);
        }
    }
}