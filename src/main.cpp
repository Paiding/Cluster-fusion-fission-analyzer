#include <fstream>
#include <vector>
#include <iostream>

#include <ctime>
#include <string>

#include "dbscan.h"
#include "analy.h"
#include "split_analy.h"

#include <unistd.h>


int main(int argc, char **argv)
{
    //clustering();

    //splitTraceBack();

    //outputDensity();

    //v6_intersecTraceBack(3300,5000);
    //v7_intersecTraceBack(0,5000);
    
    //intersecTraceBack();
    
    //splitComputeAverageAlpha(400,0,5000,1);
    //splitComputeAverageBeta(400,0,1000,1);

    //splitComputeAverageRate(400,0,5000,1);
    
    //splitCoarseGrainedAverage();

    //countSourceGain(500,1000,100);

    //outputAvgSize(2800,3500,51200,8800);

    inhomoAnaly(3);


    return 0;
}