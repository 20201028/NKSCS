#include "DataGraph.h"
#include <algorithm>
#include <memory>
DataGraph::DataGraph() : vNum(0), eNum(0) {}


void DataGraph::addEdgeNoMatinC(int src, int dst)
{
	if (src != dst)
	{
		// cout<<src<<" "<<dst<<endl;
		unordered_map<int, int>::iterator it; // unordered_map is the hashtable
		int seq1, seq2;						  // index of src/dst in the AdjList, AdjList[index] is the neighborhood of src/dst

		it = id2seq.find(src);
		if (it == id2seq.end()) // cannot find the node src, i.e., add the new node to each storage structure
		{
			seq1 = AdjList.size();
			id2seq[src] = seq1;
			seq2id[seq1] = src;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq1] = nodeAtt[src];
		}
		else
		{
			seq1 = it->second;
		}

		it = id2seq.find(dst);
		if (it == id2seq.end())
		{
			seq2 = AdjList.size();
			id2seq[dst] = seq2;
			seq2id[seq2] = dst;
			vector<int> newVertex;
			AdjList.push_back(newVertex);
			// seq2att[seq2] = nodeAtt[dst];
		}
		else
		{
			seq2 = it->second;
		}
		int flag;

		flag = 0;
		for (int i = 0; i < AdjList[seq1].size(); i++)
		{
			if (AdjList[seq1][i] == seq2)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq1].push_back(seq2); // if edge (src, dst) is not recorded

		flag = 0;
		for (int i = 0; i < AdjList[seq2].size(); i++)
		{
			if (AdjList[seq2][i] == seq1)
			{
				flag = 1;
				break;
			}
		}
		if (flag == 0)
			AdjList[seq2].push_back(seq1);
	}
}