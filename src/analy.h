#ifndef ANALY_H
#define ANALY_H

#include <fstream>
#include <vector>
#include <iostream>
#include "dbscan.h"
#include <ctime>
#include <string>
#include <algorithm>
#include <sstream>
#include "constant.h"

void readFile(std::string fileName, vector<Point>& points, int nline, vector<double>& borderList);

int arrangeCluster(vector<Point>& points, vector<int>& Size, std::vector< std::vector<int> > Index,int nPoints);

void densityAnaly(int nMonomer,int maxTypeNum, std::vector<int> clusterSize,std::vector<double>& clusterDensity);

void clustering();

void outputDensity();

void countSourceGain(int init_ser,int end_ser,int max_sections);

void outputAvgSize(int init_ser,int end_ser,int n_all,int maxTypeNum);

void inhomoAnaly(int edgeNumber);
#endif