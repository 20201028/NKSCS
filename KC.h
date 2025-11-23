#pragma once
#include "Utility.h"
#include "KTruss.h"
// #include "Basic.h"
// #include "KScore.h"
// #include "Index.h"
map<Weight, vector<vector<int>>, greater<Weight>> KC(string dataset, vector<int> posKws, 
    int k, chrono::high_resolution_clock::time_point start_time);