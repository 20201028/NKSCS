#pragma once
#include "Utility.h"
#include "KTruss.h"
#include "Basic.h"
#include "KScore.h"
// #include "Index.h"
map<Weight, vector<vector<int>>, greater<Weight>> Best_First(int method, string dataset, vector<int> posKws, vector<int> negKws, 
    int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);
    map<Weight, vector<vector<int>>, greater<Weight>> Best_First_S(int method, string dataset, vector<int> posKws, vector<int> negKws, 
    int k, Weight sthd1, Weight sthd2, chrono::high_resolution_clock::time_point start_time);