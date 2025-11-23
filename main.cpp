#include "Basic.h"
#include "Best.h"
#include "KC.h"
#include <sstream>

vector<int> parseKValues(const string &line)
{
    vector<int> k_values;
    stringstream ss(line);
    string item;
    while (getline(ss, item, ','))
    {
        k_values.push_back(stoi(item));
    }
    return k_values;
}

// void parsePkw(const std::string &line, int &pn, int &nn)
// {
//     std::stringstream ss(line);
//     std::string item;
//     std::vector<int> items;

//     while (std::getline(ss, item, ',')) {
//         items.push_back(std::stoi(item));
//     }

//     if (items.size() >= 2) {
//         pn = items[0]; // 存逗号前的数
//         nn = items[1]; // 存逗号后的数
//     } else {
//         std::cerr << "Invalid input format" << std::endl;
//     }
// }
int parsePkw(const std::string &line)
{
    std::stringstream ss(line);
    std::string item;

    getline(ss, item, ',');
    return std::stoi(item);
}

Weight parseThreshold(const string &line)
{
    stringstream ss(line);
    string numerator_str, denominator_str;
    getline(ss, numerator_str, '/');
    getline(ss, denominator_str, '/');
    int numerator = stoi(numerator_str);
    int denominator = stoi(denominator_str);
    return Weight(numerator, denominator);
}
vector<vector<int>> parseQueries(ifstream &queryFile)
{
    vector<vector<int>> queries;
    string line;

    while (getline(queryFile, line))
    {
        vector<int> keywords;
        stringstream ss(line);
        string item;
        while (getline(ss, item, ','))
        {
            keywords.push_back(stoi(item));
        }

        if (!keywords.empty())
        {
            queries.emplace_back(keywords);
        }
    }

    return queries;
}
vector<pair<vector<int>, vector<int>>> parseQueries(ifstream &queryFile, int pn)
{
    vector<pair<vector<int>, vector<int>>> queries;
    string line;

    while (getline(queryFile, line))
    {
        vector<int> keywords;
        stringstream ss(line);
        string item;
        while (getline(ss, item, ','))
        {
            keywords.push_back(stoi(item));
        }

        if (!keywords.empty())
        {
            vector<int> posKws(keywords.begin(), keywords.begin() + pn);
            vector<int> negKws(keywords.begin() + pn, keywords.end());
            queries.emplace_back(posKws, negKws);
        }
    }

    return queries;
}
vector<pair<vector<int>, vector<int>>> parseQueries(ifstream &queryFile, int pn, int nn)
{
    vector<pair<vector<int>, vector<int>>> queries;
    string line;

    while (getline(queryFile, line))
    {
        vector<int> keywords;
        stringstream ss(line);
        string item;
        while (getline(ss, item, ','))
        {
            keywords.push_back(stoi(item));
        }

        if (!keywords.empty())
        {
            vector<int> posKws(keywords.begin(), keywords.begin() + pn);
            vector<int> negKws(keywords.begin() + pn, keywords.begin() + pn + nn);
            queries.emplace_back(posKws, negKws);
        }
    }

    return queries;
}
vector<vector<int>> parseQueries2(ifstream &queryFile)
{
    vector<vector<int>> queries;
    string line;

    while (getline(queryFile, line))
    {
        // cout << line << endl;
        vector<int> keywords;
        stringstream ss(line);
        string item;
        while (getline(ss, item, ','))
        {
            keywords.push_back(stoi(item));
        }

        if (!keywords.empty())
        {
            queries.emplace_back(keywords);
        }
    }

    return queries;
}
unordered_map<int, set<int>> LoadGraph(string dataset)
{
    unordered_map<int, set<int>> datagraph;
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
        datagraph[src].insert(dst);
        datagraph[dst].insert(src);
    }

    cout << "Loaded dataset successfully!" << endl;
    return move(datagraph);
}
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

//         // id2att.emplace(ver, pair<bool, pair<Weight, vector<int>>>(true, kscoreAtt));
//         id2att.emplace(ver, att);
//     }
//     ain.close();
//     cout << "Loaded attribute successfully!" << endl;
//     return move(id2att);
// }
DataGraph LoadGraph1(string dataset)
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
        datagraph.addEdgeNoMatinC(src, dst);
    }

    cout << "Loaded dataset successfully!" << endl;
    return move(datagraph);
}
unordered_map<int, set<int>> get_subgraph(const vector<int> &nodes, const unordered_map<int, set<int>> &graph)
{
    unordered_map<int, set<int>> induced_subgraph;

    for (int u : nodes)
    {
        // 如果当前节点存在于原图中
        if (graph.find(u) != graph.end())
        {
            set<int> neighbors;

            // 遍历该节点的所有邻居，并筛选出也在子图中的节点
            for (int v : graph.at(u))
            {
                if (find(nodes.begin(), nodes.end(), v) != nodes.end())
                {
                    neighbors.insert(v);
                }
            }

            induced_subgraph[u] = move(neighbors); // 构建导出子图
        }
    }

    return induced_subgraph;
}
double get_aks(vector<int> &subgraph_nodes, const unordered_map<int, vector<int>> &att, set<int> &kw)
{
    double aks = 0.0;

    for (int i = 0; i < subgraph_nodes.size(); ++i)
    {
        int u = subgraph_nodes[i];
        const vector<int> &att_u = att.at(u);

        set<int> set_u(att_u.begin(), att_u.end());

        set<int> intersection;
        set_intersection(set_u.begin(), set_u.end(),
                         kw.begin(), kw.end(),
                         inserter(intersection, intersection.begin()));

        aks = (double)intersection.size() / (set_u.size() + kw.size() - intersection.size());
    }
    return aks;
}
double get_apk(vector<int> &subgraph_nodes, const unordered_map<int, vector<int>> &att, set<int> &kw)
{
    double apk = 0.0;
    for (int i = 0; i < subgraph_nodes.size(); ++i)
    {
        for (int j = i + 1; j < subgraph_nodes.size(); ++j)
        {
            int u = subgraph_nodes[i], v = subgraph_nodes[j];
            const vector<int> &att_u = att.at(u);
            const vector<int> &att_v = att.at(v);

            set<int> set_u(att_u.begin(), att_u.end());
            set<int> set_v(att_v.begin(), att_v.end());

            set<int> inter_uv;
            set_intersection(set_u.begin(), set_u.end(),
                             set_v.begin(), set_v.end(),
                             inserter(inter_uv, inter_uv.begin()));

            set<int> inter_q_uv;
            set_intersection(kw.begin(), kw.end(),
                             inter_uv.begin(), inter_uv.end(),
                             inserter(inter_q_uv, inter_q_uv.begin()));

            set<int> union_q_uv;
            set_union(kw.begin(), kw.end(),
                      inter_uv.begin(), inter_uv.end(),
                      inserter(union_q_uv, union_q_uv.begin()));

            apk += (double)inter_q_uv.size() / union_q_uv.size();
        }
    }
    return apk;
}
std::vector<Weight> getFractionValues(const std::string &line)
{
    std::vector<Weight> values;
    std::stringstream ss(line);
    std::string fractionStr;

    // 按逗号分割
    while (std::getline(ss, fractionStr, ','))
    {
        // 查找斜杠位置
        size_t slashPos = fractionStr.find('/');
        if (slashPos != std::string::npos)
        {
            try
            {
                int num = std::stoi(fractionStr.substr(0, slashPos));
                int den = std::stoi(fractionStr.substr(slashPos + 1));
                Weight w(num, den);
                values.push_back(w);
            }
            catch (const std::exception &e)
            {
                // 忽略解析错误
            }
        }
    }

    return values;
}
int main(int argc, char *argv[])
{

    string dataset = argv[1];
    string type = argv[2];
    int method = stoi(argv[3]);
    max_time = stof(argv[4]);
    string QueryFileName = "DataGraph/" + dataset + "/" + type + "_test.txt";

    if (type == "i")
    {
        string base_path = "DataGraph/" + dataset + "/index/";
        string StructFileName = base_path + "k.txt";
        // string StructFileName = base_path + to_string(k) + ".txt"; // src, dst, kscore
        ifstream sin(StructFileName.c_str());
        auto start_time = std::chrono::high_resolution_clock::now();

        if (!sin)
        {
            cout << "Fail to read " << StructFileName << ". Creating..." << endl;
            // truss decomposition
            // read graph

            // construct index
            // map<int, shared_ptr<ETNode>> kScore = AdvancedIndex(dg);
            Index index(dataset, base_path);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start_time;
            double search_time = duration.count();
            cout << "Construct index finish. Cost " << search_time << "s" << endl;
            // index.saveIndexToFile(StructFileName);
            // map<Weight, vector<Edge>, greater<Weight>> sim_edges = index.sim_index[k];
            // unordered_map<int, vector<int>> att_vertex = index.att_index[k];
            // unordered_map<int, vector<int>> ver_att = index.ver_index[k];
        }
        else
        {
            cout << "Exists index." << endl;
        }
        sin.close();
    }
    else if (type == "k")
    {
        ifstream queryFile(QueryFileName.c_str());
        string line;
        getline(queryFile, line);
        auto k_values = parseKValues(line);
        getline(queryFile, line);
        auto thd1 = parseThreshold(line);
        getline(queryFile, line);
        auto thd2 = parseThreshold(line);
        getline(queryFile, line);
        auto pn = parsePkw(line);
        auto queries = parseQueries(queryFile, pn);
        queryFile.close();
        // int methods = 9; // 固定方法编号或通过参数传入

        for (int k : k_values)
        {

            // auto end = std::chrono::high_resolution_clock::now();
            cout << "Running with k=" << k << endl;
            double total_time = 0.0;
            double total_time_index = 0.0;
            int query_count = 0;
            for (const auto &[posKws, negKws] : queries)
            {
                cout << "posKws: ";
                for (auto p : posKws)
                    cout << p << " ";
                cout << ", negKws: ";
                for (auto n : negKws)
                    cout << n << " ";
                cout << endl;
                auto start = std::chrono::high_resolution_clock::now();
                if (method == 0)
                {
                    Enumerate(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 1)
                {
                    AvdEnumerate(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 2)
                {
                    AvdEnumerateIncG(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 3)
                {
                    AvdEnumerate1(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 4)
                {
                    Best_First(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 5)
                {
                    AvdEnumerate(1, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_NoConntAndSize
                }
                else if (method == 6)
                {
                    AvdEnumerate(2, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_Connt
                }
                else if (method == 7)
                {
                    AvdEnumerate(3, dataset, posKws, negKws, k, thd1, thd2, start); ////AvdEnumerate_Size
                }
                else if (method == 8)
                {
                    Best_First_S(0, dataset, posKws, negKws, k, thd1, thd2, start);
                    // NoSpecial(0, dataset, posKws, negKws, k, thd1, thd2, start); // no special vertex Greedy
                    // Index_Enumerate(-1, dataset, posKws, negKws, k, thd1, thd2, start);
                    // AvdEnumerate2(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                // else if (method == 9)
                // {
                //     GE(dataset, posKws, negKws, k, thd1, thd2);
                // }
                // end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // double search_time = duration.count();
                // if (search_time > max_time)
                // {
                //     break;
                // }
                query_count++;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                cout << "search time: " << duration.count() << "s" << endl;
                double search_time = min(duration.count(), max_time);
                cout << "search time: " << search_time << "s" << endl;
                total_time += search_time;
                // auto start = std::chrono::high_resolution_clock::now();
                // vector<pair<Weight, vector<int>>> result = Enumerate(methods, dataset, posKws, negKws, k, thd1, thd2); // basic
                // auto end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // cout << "basic search time: " << duration.count() << "s" << endl;
                // auto start2 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end2 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration2 = end2 - start2;
                // cout << "bs search time: " << duration2.count() << "s" << endl;
                // auto start4 = std::chrono::high_resolution_clock::now();
                // AvdEnumerateIncG(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end4 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration4 = end4 - start4;
                // cout << "inc bs search time: " << duration4.count() << "s" << endl;
                // auto start3 = std::chrono::high_resolution_clock::now();
                // Best_First(dataset, posKws, negKws, k, thd1, thd2);
                // auto end3 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration3 = end3 - start3;
                // cout << "best search time: " << duration3.count() << "s" << endl;
                // auto start1 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate1(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end1 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration1 = end1 - start1;
                // cout << "st search time: " << duration1.count() << "s" << endl;
            }

            // std::chrono::duration<double> duration = end - start;
            // double search_time = duration.count();
            // cout << "total search time: " << search_time / query_count << "s" << endl;
            cout << "total search time: " << total_time / query_count << "s" << endl;
        }
    }
    else if (type == "e")
    {
        unordered_map<int, set<int>> graph = LoadGraph(dataset);
        unordered_map<int, vector<int>> att = LoadAtt(dataset);
        ifstream queryFile(QueryFileName.c_str());
        string line;
        getline(queryFile, line);
        auto k = stoi(line);
        getline(queryFile, line);
        auto thd1 = parseThreshold(line);
        getline(queryFile, line);
        auto thd2 = parseThreshold(line);
        getline(queryFile, line);
        auto pn = parsePkw(line);
        auto queries = parseQueries(queryFile, pn);
        queryFile.close();

        double p_aks_total = 0.0;
        double n_aks_total = 0.0;
        double apa_total = 0.0;
        double p_apk_total = 0.0;
        double n_apk_total = 0.0;
        double ed_total = 0.0;
        double nsize_total = 0.0;

        int result_num = 0;
        auto start = std::chrono::high_resolution_clock::now();
        for (const auto &[posKws, negKws] : queries)
        {
            cout << "posKws: ";
            for (auto p : posKws)
                cout << p << " ";
            cout << ", negKws: ";
            for (auto n : negKws)
                cout << n << " ";
            cout << endl;
            // vector<int> negKws;
            map<Weight, vector<vector<int>>, greater<Weight>> result;
            if (method == 1)
            {
                result = Best_First(1, dataset, posKws, negKws, k, thd1, thd2, start);
            }
            else if (method == 2)
            {
                result = Best_First(2, dataset, posKws, negKws, k, thd1, thd2, start);
            }
            else
            {
                result = KC(dataset, posKws, k, start);
            }
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            double search_time = duration.count();
            cout << "Search time: " << search_time << "s" << "max time: " << max_time << endl;
            if (search_time > max_time)
            {
                break;
            }
            if (result.empty())
            {
                continue;
            }
            for (auto &subgraph_nodes : result.begin()->second)
            {
                result_num++;

                int n = subgraph_nodes.size();
                int neg_num = 0;
                for (int u : subgraph_nodes)
                {
                    const vector<int> &att_u = att.at(u);
                    vector<int> com;
                    set_intersection(att_u.begin(), att_u.end(), negKws.begin(), negKws.end(), back_inserter(com));
                    if (com.empty())
                        continue;
                    neg_num++;
                }
                if (neg_num > 0)
                {
                    double nsize = neg_num / (double)n;
                    cout << "nsize: " << nsize << endl;
                    nsize_total += nsize;
                }
                int perfect_edge_count = n * (n - 1) / 2;
                set<int> pkw(posKws.begin(), posKws.end());
                set<int> nkw(negKws.begin(), negKws.end());
                // aks = 累加（（查询关键词 交 节点属性） /（查询关键词 并 节点属性）） / 节点数
                // double aks = 0.0;
                // set<int> kw(posKws.begin(), posKws.end());
                // for (int i = 0; i < subgraph_nodes.size(); ++i)
                // {
                //     int u = subgraph_nodes[i];
                //     const vector<int> &att_u = att.at(u);

                //     set<int> set_u(att_u.begin(), att_u.end());

                //     set<int> intersection;
                //     set_intersection(set_u.begin(), set_u.end(),
                //                      kw.begin(), kw.end(),
                //                      inserter(intersection, intersection.begin()));

                //     aks = (double)intersection.size() / (set_u.size() + kw.size() - intersection.size());
                // }
                double p_aks = get_aks(subgraph_nodes, att, pkw);
                if (p_aks > 0)
                {
                    p_aks /= subgraph_nodes.size();
                    cout << "positive aks: " << p_aks << endl;
                    p_aks_total += p_aks;
                }
                double n_aks = get_aks(subgraph_nodes, att, nkw);
                if (n_aks > 0)
                {
                    n_aks /= subgraph_nodes.size();
                    cout << "negative aks: " << n_aks << endl;
                    n_aks_total += n_aks;
                }
                // 平均成对属性相似性apa = 累加任意两节点(（u节点属性 交 v节点属性） /（u节点属性 并 v节点属性）) / 完美边数
                // vector<int> nodes(subgraph_nodes.begin(), subgraph_nodes.end());
                double apa = 0.0;
                for (int i = 0; i < subgraph_nodes.size(); ++i)
                {
                    for (int j = i + 1; j < subgraph_nodes.size(); ++j)
                    {
                        int u = subgraph_nodes[i], v = subgraph_nodes[j];
                        const vector<int> &att_u = att.at(u);
                        const vector<int> &att_v = att.at(v);

                        set<int> set_u(att_u.begin(), att_u.end());
                        set<int> set_v(att_v.begin(), att_v.end());

                        set<int> intersection;
                        set_intersection(set_u.begin(), set_u.end(),
                                         set_v.begin(), set_v.end(),
                                         inserter(intersection, intersection.begin()));

                        apa += (double)intersection.size() / (set_u.size() + set_v.size() - intersection.size());
                    }
                }
                if (apa > 0)
                {
                    apa /= perfect_edge_count;
                    cout << "apa: " << apa << endl;
                    apa_total += apa;
                }
                // 平均成对关键词相似性apk = 累加任意两节点(（查询关键词 交(u节点属性 交 v节点属性)）/（查询关键词 并(u节点属性 交 v节点属性)）) / 边数
                // double apk = 0.0;
                // for (int i = 0; i < subgraph_nodes.size(); ++i)
                // {
                //     for (int j = i + 1; j < subgraph_nodes.size(); ++j)
                //     {
                //         int u = subgraph_nodes[i], v = subgraph_nodes[j];
                //         const vector<int> &att_u = att.at(u);
                //         const vector<int> &att_v = att.at(v);

                //         set<int> set_u(att_u.begin(), att_u.end());
                //         set<int> set_v(att_v.begin(), att_v.end());

                //         set<int> inter_uv;
                //         set_intersection(set_u.begin(), set_u.end(),
                //                          set_v.begin(), set_v.end(),
                //                          inserter(inter_uv, inter_uv.begin()));

                //         set<int> inter_q_uv;
                //         set_intersection(kw.begin(), kw.end(),
                //                          inter_uv.begin(), inter_uv.end(),
                //                          inserter(inter_q_uv, inter_q_uv.begin()));

                //         set<int> union_q_uv;
                //         set_union(kw.begin(), kw.end(),
                //                   inter_uv.begin(), inter_uv.end(),
                //                   inserter(union_q_uv, union_q_uv.begin()));

                //         apk += (double)inter_q_uv.size() / union_q_uv.size();
                //     }
                // }
                double p_apk = get_apk(subgraph_nodes, att, pkw);
                if (p_apk > 0)
                {
                    p_apk /= perfect_edge_count;
                    cout << "positive apk: " << p_apk << endl;
                    p_apk_total += p_apk;
                }
                double n_apk = get_apk(subgraph_nodes, att, nkw);
                if (n_apk > 0)
                {
                    n_apk /= perfect_edge_count;
                    cout << "negative apk: " << n_apk << endl;
                    n_apk_total += n_apk;
                }
                // ed = 边数 / 完美边数
                unordered_map<int, set<int>> induced_subgraph = get_subgraph(subgraph_nodes, graph);
                int edge_count = 0;
                for (const auto &[u, neighbors] : induced_subgraph)
                {
                    for (int v : neighbors)
                    {
                        if (u < v)
                            edge_count++; // 避免重复计数无向边
                    }
                }
                double ed = (double)edge_count / perfect_edge_count;
                cout << "ed: " << ed << endl;
                ed_total += ed;
            }
            // cout << "search time: " << search_time << "s" << endl;
        }
        cout << "total positive aks: " << p_aks_total / (double)result_num << endl;
        cout << "total negative aks: " << n_aks_total / (double)result_num << endl;
        cout << "total apa: " << apa_total / (double)result_num << endl;
        cout << "total positive apk: " << p_apk_total / (double)result_num << endl;
        cout << "total negative apk: " << n_apk_total / (double)result_num << endl;
        cout << "total ed: " << ed_total / (double)result_num << endl;
        cout << "total nsize: " << nsize_total / (double)result_num << endl;
    }
    else if (type == "s")
    {
        ifstream queryFile(QueryFileName.c_str());
        string line;
        getline(queryFile, line);
        auto k = stoi(line);
        getline(queryFile, line);
        vector<Weight> thds1 = getFractionValues(line);
        getline(queryFile, line);
        vector<Weight> thds2 = getFractionValues(line);
        getline(queryFile, line);
        auto pn = parsePkw(line);
        auto queries = parseQueries(queryFile, pn);
        queryFile.close();
        // int methods = 9; // 固定方法编号或通过参数传入

        for (int i = 0; i < thds1.size(); i++)
        {
            Weight thd1 = thds1[i];
            Weight thd2 = thds2[i];
            // auto end = std::chrono::high_resolution_clock::now();
            cout << "Running with size = [" << thd1.numerator << '/' << thd1.denominator << ',' << thd2.numerator << '/' << thd2.denominator << ']' << endl;
            double total_time = 0.0;
            double total_time_index = 0.0;
            int query_count = 0;
            for (const auto &[posKws, negKws] : queries)
            {
                cout << "posKws: ";
                for (auto p : posKws)
                    cout << p << " ";
                cout << ", negKws: ";
                for (auto n : negKws)
                    cout << n << " ";
                cout << endl;
                auto start = std::chrono::high_resolution_clock::now();
                if (method == 0)
                {
                    Enumerate(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 1)
                {
                    AvdEnumerate(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 2)
                {
                    AvdEnumerateIncG(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 3)
                {
                    AvdEnumerate1(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 4)
                {
                    Best_First(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 5)
                {
                    AvdEnumerate(1, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_NoConntAndSize
                }
                else if (method == 6)
                {
                    AvdEnumerate(2, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_Connt
                }
                else if (method == 7)
                {
                    AvdEnumerate(3, dataset, posKws, negKws, k, thd1, thd2, start); ////AvdEnumerate_Size
                }
                else if (method == 8)
                {
                    Best_First_S(0, dataset, posKws, negKws, k, thd1, thd2, start);
                    // NoSpecial(0, dataset, posKws, negKws, k, thd1, thd2, start); // no special vertex Greedy
                    // Index_Enumerate(-1, dataset, posKws, negKws, k, thd1, thd2, start);
                    // AvdEnumerate2(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                // else if (method == 9)
                // {
                //     GE(dataset, posKws, negKws, k, thd1, thd2);
                // }
                // end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // double search_time = duration.count();
                // if (search_time > max_time)
                // {
                //     break;
                // }
                query_count++;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                cout << "search time: " << duration.count() << "s" << endl;
                double search_time = min(duration.count(), max_time);
                cout << "search time: " << search_time << "s" << endl;
                total_time += search_time;
                // auto start = std::chrono::high_resolution_clock::now();
                // vector<pair<Weight, vector<int>>> result = Enumerate(methods, dataset, posKws, negKws, k, thd1, thd2); // basic
                // auto end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // cout << "basic search time: " << duration.count() << "s" << endl;
                // auto start2 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end2 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration2 = end2 - start2;
                // cout << "bs search time: " << duration2.count() << "s" << endl;
                // auto start4 = std::chrono::high_resolution_clock::now();
                // AvdEnumerateIncG(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end4 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration4 = end4 - start4;
                // cout << "inc bs search time: " << duration4.count() << "s" << endl;
                // auto start3 = std::chrono::high_resolution_clock::now();
                // Best_First(dataset, posKws, negKws, k, thd1, thd2);
                // auto end3 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration3 = end3 - start3;
                // cout << "best search time: " << duration3.count() << "s" << endl;
                // auto start1 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate1(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end1 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration1 = end1 - start1;
                // cout << "st search time: " << duration1.count() << "s" << endl;
            }

            // std::chrono::duration<double> duration = end - start;
            // double search_time = duration.count();
            // cout << "total search time: " << search_time / query_count << "s" << endl;
            cout << "total search time: " << total_time / query_count << "s" << endl;
        }
    }
    else if (type == "p")
    {
        ifstream queryFile(QueryFileName.c_str());
        string line;
        getline(queryFile, line);
        auto k = stoi(line);
        getline(queryFile, line);
        auto thd1 = parseThreshold(line);
        getline(queryFile, line);
        auto thd2 = parseThreshold(line);
        getline(queryFile, line);
        auto pns = parseKValues(line);
        auto queries = parseQueries(queryFile);

        // int methods = 9; // 固定方法编号或通过参数传入

        for (int i = 0; i < pns.size(); i++)
        {
            int pn = pns[i];

            // auto end = std::chrono::high_resolution_clock::now();
            cout << "Running with pn =" << pn << endl;
            double total_time = 0.0;
            double total_time_index = 0.0;
            int query_count = 0;
            map<string, double> query_str;
            for (auto &kws : queries)
            {
                vector<int> posKws;
                vector<int> negKws;
                for (int i = 0; i < pn; i++)
                {
                    posKws.push_back(kws[i]);
                }
                cout << "posKws: ";
                for (auto p : posKws)
                    cout << p << " ";
                cout << ", negKws: ";
                for (auto n : negKws)
                    cout << n << " ";
                cout << endl;

                string s;
                if (posKws.size() == 1)
                {
                    s += to_string(posKws[posKws.size() - 1]);
                }
                else if (posKws.size() > 1)
                {
                    for (int i = 0; i < posKws.size() - 1; i++)
                    {
                        s += to_string(posKws[i]) + ",";
                    }
                    s += to_string(posKws[posKws.size() - 1]);
                }

                if (query_str.count(s))
                {
                    cout << "Repeat Query" << endl;
                    total_time += query_str[s];
                    query_count++;
                    continue;
                }
                auto start = std::chrono::high_resolution_clock::now();
                if (method == 0)
                {
                    Enumerate(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 1)
                {
                    AvdEnumerate(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 2)
                {
                    AvdEnumerateIncG(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 3)
                {
                    AvdEnumerate1(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 4)
                {
                    Best_First(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 5)
                {
                    AvdEnumerate(1, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_NoConntAndSize
                }
                else if (method == 6)
                {
                    AvdEnumerate(2, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_Connt
                }
                else if (method == 7)
                {
                    AvdEnumerate(3, dataset, posKws, negKws, k, thd1, thd2, start); ////AvdEnumerate_Size
                }
                else if (method == 8)
                {
                    Best_First_S(0, dataset, posKws, negKws, k, thd1, thd2, start);
                    // NoSpecial(0, dataset, posKws, negKws, k, thd1, thd2, start); // no special vertex Greedy
                    // Index_Enumerate(-1, dataset, posKws, negKws, k, thd1, thd2, start);
                    // AvdEnumerate2(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                // else if (method == 9)
                // {
                //     GE(dataset, posKws, negKws, k, thd1, thd2);
                // }
                // end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // double search_time = duration.count();
                // if (search_time > max_time)
                // {
                //     break;
                // }
                query_count++;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                cout << "search time: " << duration.count() << "s" << endl;
                double search_time = min(duration.count(), max_time);
                query_str[s] = search_time;
                cout << "search time: " << search_time << "s" << endl;
                total_time += search_time;
                // auto start = std::chrono::high_resolution_clock::now();
                // vector<pair<Weight, vector<int>>> result = Enumerate(methods, dataset, posKws, negKws, k, thd1, thd2); // basic
                // auto end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // cout << "basic search time: " << duration.count() << "s" << endl;
                // auto start2 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end2 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration2 = end2 - start2;
                // cout << "bs search time: " << duration2.count() << "s" << endl;
                // auto start4 = std::chrono::high_resolution_clock::now();
                // AvdEnumerateIncG(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end4 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration4 = end4 - start4;
                // cout << "inc bs search time: " << duration4.count() << "s" << endl;
                // auto start3 = std::chrono::high_resolution_clock::now();
                // Best_First(dataset, posKws, negKws, k, thd1, thd2);
                // auto end3 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration3 = end3 - start3;
                // cout << "best search time: " << duration3.count() << "s" << endl;
                // auto start1 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate1(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end1 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration1 = end1 - start1;
                // cout << "st search time: " << duration1.count() << "s" << endl;
            }

            // std::chrono::duration<double> duration = end - start;
            // double search_time = duration.count();
            // cout << "total search time: " << search_time / query_count << "s" << endl;
            cout << "total search time: " << total_time / query_count << "s" << endl;
        }
        queryFile.close();
    }
    else if (type == "n")
    {
        ifstream queryFile(QueryFileName.c_str());
        string line;
        getline(queryFile, line);
        auto k = stoi(line);
        getline(queryFile, line);
        auto thd1 = parseThreshold(line);
        getline(queryFile, line);
        auto thd2 = parseThreshold(line);
        getline(queryFile, line);
        auto pn = stoi(line);
        getline(queryFile, line);
        auto nns = parseKValues(line);
        auto queries = parseQueries(queryFile);

        // int methods = 9; // 固定方法编号或通过参数传入

        for (int i = 0; i < nns.size(); i++)
        {
            int nn = nns[i];
            // auto queries = parseQueries(queryFile, pn, nn);
            // auto end = std::chrono::high_resolution_clock::now();
            cout << "Running with nn =" << nn << endl;
            double total_time = 0.0;
            double total_time_index = 0.0;
            int query_count = 0;
            map<string, double> query_str;
            for (auto &kws : queries)
            {
                vector<int> posKws;
                vector<int> negKws;
                for (int i = 0; i < pn; i++)
                {
                    posKws.push_back(kws[i]);
                }
                for (int i = pn; i < pn + nn; i++)
                {
                    negKws.push_back(kws[i]);
                }
                cout << "posKws: ";
                for (auto p : posKws)
                    cout << p << " ";
                cout << ", negKws: ";
                for (auto n : negKws)
                    cout << n << " ";
                cout << endl;
                string s = to_string(posKws[0]);
                if (negKws.size() != 0)
                {
                    for (int i = 0; i < negKws.size(); i++)
                    {
                        s += "," + to_string(negKws[i]);
                    }
                }

                if (query_str.count(s))
                {
                    cout << "Repeat Query" << endl;
                    total_time += query_str[s];
                    query_count++;
                    continue;
                }
                auto start = std::chrono::high_resolution_clock::now();
                if (method == 0)
                {
                    Enumerate(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 1)
                {
                    AvdEnumerate(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 2)
                {
                    AvdEnumerateIncG(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 3)
                {
                    AvdEnumerate1(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 4)
                {
                    Best_First(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 5)
                {
                    AvdEnumerate(1, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_NoConntAndSize
                }
                else if (method == 6)
                {
                    AvdEnumerate(2, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_Connt
                }
                else if (method == 7)
                {
                    AvdEnumerate(3, dataset, posKws, negKws, k, thd1, thd2, start); ////AvdEnumerate_Size
                }
                else if (method == 8)
                {
                    Best_First_S(0, dataset, posKws, negKws, k, thd1, thd2, start);
                    // NoSpecial(0, dataset, posKws, negKws, k, thd1, thd2, start); // no special vertex Greedy
                    // Index_Enumerate(-1, dataset, posKws, negKws, k, thd1, thd2, start);
                    // AvdEnumerate2(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                // else if (method == 9)
                // {
                //     GE(dataset, posKws, negKws, k, thd1, thd2);
                // }
                // end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // double search_time = duration.count();
                // if (search_time > max_time)
                // {
                //     break;
                // }
                query_count++;
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                cout << "search time: " << duration.count() << "s" << endl;
                double search_time = min(duration.count(), max_time);
                query_str[s] = search_time;
                cout << "search time: " << search_time << "s" << endl;
                total_time += search_time;
                // auto start = std::chrono::high_resolution_clock::now();
                // vector<pair<Weight, vector<int>>> result = Enumerate(methods, dataset, posKws, negKws, k, thd1, thd2); // basic
                // auto end = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration = end - start;
                // cout << "basic search time: " << duration.count() << "s" << endl;
                // auto start2 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end2 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration2 = end2 - start2;
                // cout << "bs search time: " << duration2.count() << "s" << endl;
                // auto start4 = std::chrono::high_resolution_clock::now();
                // AvdEnumerateIncG(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end4 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration4 = end4 - start4;
                // cout << "inc bs search time: " << duration4.count() << "s" << endl;
                // auto start3 = std::chrono::high_resolution_clock::now();
                // Best_First(dataset, posKws, negKws, k, thd1, thd2);
                // auto end3 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration3 = end3 - start3;
                // cout << "best search time: " << duration3.count() << "s" << endl;
                // auto start1 = std::chrono::high_resolution_clock::now();
                // AvdEnumerate1(methods, dataset, posKws, negKws, k, thd1, thd2);
                // auto end1 = std::chrono::high_resolution_clock::now();
                // std::chrono::duration<double> duration1 = end1 - start1;
                // cout << "st search time: " << duration1.count() << "s" << endl;
            }

            // std::chrono::duration<double> duration = end - start;
            // double search_time = duration.count();
            // cout << "total search time: " << search_time / query_count << "s" << endl;
            cout << "total search time: " << total_time / query_count << "s" << endl;
        }
        queryFile.close();
    }
    else if (type == "d")
    {
        // dataset statistics
        DataGraph dg = LoadGraph1(dataset);
        double avg_deg = 0;
        for (auto &nei : dg.AdjList)
        {
            avg_deg += nei.size();
        }
        cout << "vNum: " << dg.AdjList.size() << endl;
        cout << "eNum: " << avg_deg / 2 << endl;
        cout << "avg_deg: " << avg_deg / dg.AdjList.size() << endl;
        auto start = std::chrono::high_resolution_clock::now();
        TrussDecomposition *td = new TrussDecomposition(dg);
        cout << "kMax: " << td->kMax << endl;
        unordered_map<int, vector<int>> id2att = LoadAtt(dataset);
        double avg_att = 0;
        unordered_set<int> att_set;
        for (auto &ia : id2att)
        {
            avg_att += ia.second.size();
            att_set.insert(ia.second.begin(), ia.second.end());
        }
        cout << "avg_att: " << avg_att / id2att.size() << endl;
        cout << "att_num: " << att_set.size() << endl;
    }
    else if (type == "v")
    {
        double temp_max_time = max_time;
        // max_time = 1000000000;
        ifstream queryFile(QueryFileName.c_str());
        string line;
        getline(queryFile, line);
        auto k = stoi(line);
        getline(queryFile, line);
        auto thd1 = parseThreshold(line);
        getline(queryFile, line);
        auto thd2 = parseThreshold(line);
        getline(queryFile, line);
        auto pn = parsePkw(line);
        auto queries = parseQueries(queryFile, pn);
        queryFile.close();
        int query_count = 0;
        double total_time = 0.0;
        auto start = std::chrono::high_resolution_clock::now();
        for (const auto &[posKws, negKws] : queries)
        {
            cout << "posKws: ";
            for (auto p : posKws)
                cout << p << " ";
            cout << ", negKws: ";
            for (auto n : negKws)
                cout << n << " ";
            cout << endl;
            // vector<int> negKws;
            map<Weight, vector<vector<int>>, greater<Weight>> result;

            if (method == 0)
            {
                result = AvdEnumerate1(9, dataset, posKws, negKws, k, thd1, thd2, start);
            }
            else if (method == 1)
            {
                result = Best_First_S(0, dataset, posKws, negKws, k, thd1, thd2, start);
            }
            // query_count++;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = end - start;
            double search_time = duration.count();
            cout << "Search time: " << search_time << "s" << "max time: " << temp_max_time << endl;

            if (search_time > temp_max_time)
            {
                // total_time = max_time;
                break;
            }
            query_count++;
            total_time = search_time;

            // cout << "search time: " << search_time << "s" << endl;
        }
        cout << "total search time: " << total_time / query_count << "s" << endl;
    }
    else if (type == "c")
    {
        // double temp_max_time = max_time;
        max_time = 1000000000;
        unordered_map<int, set<int>> graph = LoadGraph(dataset);
        ifstream queryFile(QueryFileName.c_str());
        string line;
        getline(queryFile, line);
        auto k = stoi(line);
        getline(queryFile, line);
        vector<Weight> thds1 = getFractionValues(line);
        getline(queryFile, line);
        vector<Weight> thds2 = getFractionValues(line);
        getline(queryFile, line);
        auto pn = parsePkw(line);
        auto queries = parseQueries(queryFile, pn);
        queryFile.close();
        // int methods = 9; // 固定方法编号或通过参数传入

        for (int i = 0; i < thds1.size(); i++)
        {
            Weight thd1 = thds1[i];
            Weight thd2 = thds2[i];
            // auto end = std::chrono::high_resolution_clock::now();
            cout << "Running with size = [" << thd1.numerator << '/' << thd1.denominator << ',' << thd2.numerator << '/' << thd2.denominator << ']' << endl;
            for (const auto &[posKws, negKws] : queries)
            {
                cout << "posKws: ";
                for (auto p : posKws)
                    cout << p << " ";
                cout << ", negKws: ";
                for (auto n : negKws)
                    cout << n << " ";
                cout << endl;
                map<Weight, vector<vector<int>>, greater<Weight>> result;
                auto start = std::chrono::high_resolution_clock::now();
                if (method == 0)
                {
                    Enumerate(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 1)
                {
                    AvdEnumerate(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 2)
                {
                    result = Best_First(0, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 3)
                {
                    result = AvdEnumerate1(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 4)
                {
                    result = Best_First(2, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                else if (method == 5)
                {
                    AvdEnumerate(1, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_NoConntAndSize
                }
                else if (method == 6)
                {
                    AvdEnumerate(2, dataset, posKws, negKws, k, thd1, thd2, start); // AvdEnumerate_Connt
                }
                else if (method == 7)
                {
                    AvdEnumerate(3, dataset, posKws, negKws, k, thd1, thd2, start); ////AvdEnumerate_Size
                }
                else if (method == 8)
                {
                    Best_First_S(0, dataset, posKws, negKws, k, thd1, thd2, start);
                    // NoSpecial(0, dataset, posKws, negKws, k, thd1, thd2, start); // no special vertex Greedy
                    // Index_Enumerate(-1, dataset, posKws, negKws, k, thd1, thd2, start);
                    // AvdEnumerate2(9, dataset, posKws, negKws, k, thd1, thd2, start);
                }
                std::cout << "edges = [\n";
                for (auto &subgraph_nodes : result.begin()->second)
                {
                    
                    unordered_map<int, set<int>> induced_subgraph = get_subgraph(subgraph_nodes, graph);
                    int edge_count = 0;
                    for (const auto &[u, neighbors] : induced_subgraph)
                    {
                        for (int v : neighbors)
                        {
                            if(u > v){
                                continue;
                            }
                            std::cout << "    (" << u << ", " << v << "),\n";
                        }
                    }
                    
                }
                std::cout << "]\n";

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;
                cout << "search time: " << duration.count() << "s" << endl;
            }
        }
    }

    return 0;
}