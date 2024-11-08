#ifndef SPLIT_ANALY_H
#define SPLIT_ANALY_H

#include <fstream>
#include <vector>
#include <iostream>
#include "dbscan.h"
#include <ctime>
#include <string>
#include <algorithm>
#include <sstream>
#include "constant.h"

void intersecTraceBack();

void v7_intersecTraceBack(int init_ser,int end_ser);

void v6_intersecTraceBack(int init_ser,int end_ser);

void splitComputeAverageRate(int maxTypeNum,int beginFileSer,int endFileSer,int delta_t);

void splitComputeAverageAlpha(int maxTypeNum,int beginFileSer,int endFileSer,int delta_t);

void splitComputeAverageBeta(int maxTypeNum,int beginFileSer,int endFileSer,int delta_t);

void splitCoarseGrainedAverage();

void splitCoarseGrainedAverage_withName(std::string inFilename,std::string outFilename);
#endif