#include "analy.h"

void readFile(std::string fileName, vector<Point>& points, int nline, vector<double>& borderList)
{
    std::ifstream file;
    std::string line;
    Point *p = (Point *)calloc(nline, sizeof(Point));
    file.open(fileName);
    getline(file,line);
    getline(file,line);
    for(int i = 0; i < nline; ++i)
    {
        getline(file, line);
        p[i].clusterID = UNCLASSIFIED;
        p[i].x = stof(line.substr(22, 6));
        p[i].y = stof(line.substr(30, 6));
        p[i].z = stof(line.substr(38, 6));
        points.push_back(p[i]);
    }
    free(p);
    //added for borderlist
    getline(file,line);
    std::istringstream iss(line);
    double number;
    while (iss >> number)
    {
        borderList.push_back(number);
    }
    //end add
    file.close();
}

int arrangeCluster(std::vector<Point>& points, std::vector<int>& Size, std::vector< std::vector<int> > Index,int nPoints)
{
    int nMonomer = 0;
    for (int i = 0; i < nPoints; i++)
    {
        int current_index = points[i].clusterID;
        if (current_index < 0)
        {
            nMonomer ++;
        }
        else
        {
            Size[current_index] = Size[current_index] + 1;
            Index[current_index].push_back(i);
        }
    }
    return nMonomer;
}

void densityAnaly(int nMonomer,int maxTypeNum, std::vector<int> clusterSize,std::vector<double>& clusterDensity)
{
    clusterDensity.clear();
    clusterDensity.resize(maxTypeNum,0.0);
    for (int i = 0; i < clusterSize.size(); i++)
    {
        int current_size = clusterSize[i];
        if (current_size < maxTypeNum)
        {
            clusterDensity[current_size] += DENSCOEF;
        }
        
    }
    clusterDensity[1] = nMonomer * DENSCOEF;
}

void clustering()
{
    clock_t start, finish;
    start = clock(); 
    int Natom = 32;
    int Nmole = 400;
    int nline = Natom * Nmole;
    //int nline = 3200;
    float eps = 0.55;
    int minSize = Natom + 1;
    int count = 0; 
    int current_index = 0;
    int init_fileNameIndex = 0;
    int end_fileNameIndex = 2500;
    
    std::vector<int> iSize_t0(nline, 1);
    std::vector<int> iSize_t1(nline, 1);
    std::vector<int> iNdx_t0(nline, 0);
    std::vector<int> iNdx_t1(nline, 0);
    
    std::string logFilename = std::to_string(init_fileNameIndex) + "_log.txt"; 
    std::ofstream logFileStream;
    logFileStream.open(logFilename);

    for (int i_filename = init_fileNameIndex; i_filename < end_fileNameIndex; i_filename++)
    {
    /*Read file*/
        std::string filename = "conf/conf" + std::to_string(i_filename) + ".gro";
        std::string clustIDFilename = "Index/conf" + std::to_string(i_filename) + "_id.dat";
        std::string densityFilename = "test/conf" + std::to_string(i_filename) + "_c.dat";
        
        std::vector<double> borderXYZ;
        std::vector<Point> points;
        readFile(filename,points,nline,borderXYZ);

    /*Initialize interstep vectors*/
        for (int i = 0; i < nline; i++)
        {
            iSize_t0[i] = iSize_t1[i];
            iSize_t1[i] = 1;
            iNdx_t0[i] = iNdx_t1[i];
            iNdx_t1[i] = 0; 
        }
        
    /*Clustering*/
        DBSCAN ds(minSize,eps,points,borderXYZ);
        int nclust = ds.run();
        int nMonomer = 0;
        int maxclust = 0;

    /*Analize*/
        float averageSize = 0;
        std::vector<int> deltaValue(nline, 0);
        std::cout<<"Nclust= "<<nclust<<", ";
        std::vector<int> clusterSize(nclust, 0);
        std::vector<std::vector<int>> clusterIndex(nclust,std::vector<int>());
        std::vector<double> clusterDensity;
        
        //sort the point array and arrange them to clustIndex
        /*for multiatom molecules*/
        for (int i = 0; i < Nmole; i++)
        {
            int first_atom_index = i * Natom;
            int current_index = ds.m_points[first_atom_index].clusterID;
            if (current_index < 0)
            {
                nMonomer ++;
            }
            else
            {
                clusterSize[current_index] = clusterSize[current_index] + 1;
                clusterIndex[current_index].push_back(i);
            }
        }
        
        //calculate and log basic clustering parameters for this step
        
        std::cout<<"nMono= "<<nMonomer<<", ";
        for (int i = 0; i < clusterSize.size(); i++)
        {
            averageSize += clusterSize[i];
            maxclust = (maxclust>clusterSize[i])?maxclust:(clusterSize[i]);
        }
        averageSize = averageSize / clusterSize.size();
        std::cout<<"AvgSize= "<<averageSize<<", ";
        std::cout<<"maxClustS = "<<maxclust<<"\n";
        
        logFileStream<<" "<<nclust<<" "<<nMonomer<<" "<<averageSize<<" "<<maxclust<<"\n";
        
        
        //calculate and log intermediate variables
        //densityAnaly(nMonomer,maxTypeNum,clusterSize,clusterDensity); 
        /*
        std::ofstream outStream;
        outStream.open(densityFilename);
        for (int i = 0; i < maxTypeNum; i++)
        {
            outStream<<clusterDensity[i]<<endl;
        }*/
        
        std::ofstream clustIndexFile;
        clustIndexFile.open(clustIDFilename);
        for (int i = 0; i < clusterIndex.size(); i++)
        {
            for (int j = 0; j < clusterIndex[i].size(); j++)
            {
                current_index = clusterIndex[i][j];
                clustIndexFile<<current_index;
                clustIndexFile<<",";
            }
            clustIndexFile<<"\n";
        }
        clustIndexFile.close();
    }
    logFileStream.close();
    finish = clock(); 	
    printf( "Clustering cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
}


void outputDensity()
{
    clock_t start, finish;
    start = clock(); 
    int init_ser = 500;
    int end_ser = 501;
    //int n_all = 800;
    int n_all = 400;
    //int maxTypeNum = 40;
    int maxTypeNum = 400;
    double inverseV = 6.714e-4;
    double tau = 0.02;
    double eps = 1e-9;
    int delta_t = 1;
    int current_index = 0;
    int fine_intervals = 400;
    int cg_gap = 0;
    int cg_capacity = 87;
    for (size_t fileNameSerial = init_ser; fileNameSerial < end_ser; fileNameSerial++)
    {
        std::vector<std::vector<int>> clusterIndex_t0;
        std::vector<std::vector<int>> clusterIndex_t1;
        std::vector<double> alpha(maxTypeNum,0.0);
        std::vector<std::vector<double>> beta(maxTypeNum,std::vector<double>(maxTypeNum,0.0));
        std::vector<int> iSize_t0(n_all, 1);
        std::vector<int> iSize_t1(n_all, 1);
        std::vector<int> iNdx_t0(n_all, 0);
        std::vector<int> iNdx_t1(n_all, 0);
        std::vector<double> clusterDensity(maxTypeNum,0.0);
        std::string logFilename = "frame" + std::to_string(fileNameSerial) + "Density.txt";

    /*Read File*/
        std::string clustIDFilename0 = "Index/conf" + std::to_string(fileNameSerial) + "_id.dat";

        std::ifstream file;
        std::string str;
        file.open(clustIDFilename0);
        while (std::getline(file, str)) {
            std::vector<int> row;
            std::stringstream ss(str);
            int value;
            while (ss >> value) {
                row.push_back(value);
                if (ss.peek() == ',') ss.ignore();
            }
            clusterIndex_t0.push_back(row);
        }
        file.close();
    /*Rearrangement*/
        int nclust0 = clusterIndex_t0.size();
        int nMonomer0 = n_all;//minus sizei until only monomer remain
        int nMonomer1 = n_all;
        for (int i = 0; i < nclust0; i++)
        {
            int current_size = clusterIndex_t0[i].size();
            nMonomer0 -= current_size;
            for (int j = 0; j < current_size; j++)
            {
                current_index = clusterIndex_t0[i][j];
                iSize_t0[current_index] = clusterIndex_t0[i].size();
                iNdx_t0[current_index] = i;
            }
        }
        
    /*Compute Density*/
        for (int i = 0; i < clusterIndex_t0.size(); i++)
        {
            int current_size = clusterIndex_t0[i].size();
            if (current_size < maxTypeNum)
            {
                clusterDensity[current_size] += inverseV;
            }
        }
        clusterDensity[1] = nMonomer0 * inverseV;

    /*compute cg density*/
        std::vector<double> cg_density(fine_intervals + cg_capacity,0.0);
        for (int i = 0; i < fine_intervals; i++)
        {
            cg_density[i] = clusterDensity[i];
        }
        for (int i = fine_intervals; i < fine_intervals+cg_capacity; i++)
        {
            for (int j = 0; j < cg_gap; j++)
            {
                int crindex = fine_intervals + (i-fine_intervals) * cg_gap + j;
                cg_density[i] += clusterDensity[crindex];
            }
        }
        

        std::ofstream densityOutput;
        densityOutput.open(logFilename);
        for (int i = 1; i < cg_density.size(); i++)
        {
            densityOutput<<cg_density[i]<<std::endl;
        }
        
    }
}

void countSourceGain(int init_ser,int end_ser,int max_sections)
{//count the number of source gain and end loss
    clock_t start, finish;
    start = clock(); 
    //int n_all = 800;
    int n_all = 51200;
    //int maxTypeNum = 300;
    int maxTypeNum = 4400;
    double inverseV = 1.5625e-5;
    double tau = 0.2;
    double eps = 1e-12;

    for (size_t fileNameSerial = init_ser; fileNameSerial < end_ser; fileNameSerial++)
    {
        std::vector<std::vector<int>> clusterIndex_t0;
        std::vector<std::vector<int>> clusterIndex_t1;
        std::vector<std::vector<double>> alpha(maxTypeNum,std::vector<double>(maxTypeNum,0.0));
        std::vector<std::vector<double>> beta(maxTypeNum,std::vector<double>(maxTypeNum,0.0));
        std::vector<int> iSize_t0(n_all, 1);
        std::vector<int> iSize_t1(n_all, 1);
        std::vector<int> iNdx_t0(n_all, 0);
        std::vector<int> iNdx_t1(n_all, 0);
        std::vector<double> clusterDensity(maxTypeNum,0.0);
    /*Read File*/
        std::string clustIDFilename0 = "Index/conf" + std::to_string(fileNameSerial) + "_id.dat";
        std::string clustIDFilename1 = "Index/conf" + std::to_string(fileNameSerial+1) + "_id.dat";
        //std::string clustIDFilename0 = "conf_id/conf" + std::to_string(fileNameSerial) + "_id.dat";
        //std::string clustIDFilename1 = "conf_id/conf" + std::to_string(fileNameSerial+1) + "_id.dat";
        std::ifstream file;
        std::string str;
        file.open(clustIDFilename0);
        while (std::getline(file, str)) {
            std::vector<int> row;
            std::stringstream ss(str);
            int value;
            while (ss >> value) {
                row.push_back(value);
                if (ss.peek() == ',') ss.ignore();
            }
            clusterIndex_t0.push_back(row);
        }
        file.close();
        file.open(clustIDFilename1);
        while (std::getline(file, str)) {
            std::vector<int> row;
            std::stringstream ss(str);
            int value;
            while (ss >> value) {
                row.push_back(value);
                if (ss.peek() == ',') ss.ignore();
            }
            clusterIndex_t1.push_back(row);
        }
        file.close();
    /*Rearrangement*/
        int nclust0 = clusterIndex_t0.size();
        int nclust1 = clusterIndex_t1.size();
        int nMonomer0 = n_all;//minus sizei until only monomer remain
        int nMonomer1 = n_all;
        for (int i = 0; i < nclust0; i++)
        {
            int current_size = clusterIndex_t0[i].size();
            nMonomer0 -= current_size;
            for (int j = 0; j < current_size; j++)
            {
                int current_index = clusterIndex_t0[i][j];
                iSize_t0[current_index] = clusterIndex_t0[i].size();
                iNdx_t0[current_index] = i;
            }
        }
        for (int i = 0; i < nclust1; i++)
        {
            int current_size = clusterIndex_t1[i].size();
            nMonomer1 -= current_size;
            for (int j = 0; j < current_size; j++) 
            {
                int current_index = clusterIndex_t1[i][j];
                iSize_t1[current_index] = clusterIndex_t1[i].size();
                iNdx_t1[current_index] = i;
            }
        }
    /*Compute Density*/
        for (int i = 0; i < clusterIndex_t0.size(); i++)
        {
            int current_size = clusterIndex_t0[i].size();
            if (current_size < maxTypeNum)
            {
                clusterDensity[current_size] += inverseV;
            }
        }
        clusterDensity[1] = nMonomer0 * inverseV;

    /*extra rearrangement only for intersectraceback*/
        std::vector<int> cSize_t0;
        std::vector<int> cSize_t1;
        for (int i = 0; i < nclust0; i++)
        {
            cSize_t0.push_back(clusterIndex_t0[i].size());
        }
        for (int i = 0; i < nclust1; i++)
        {
            cSize_t1.push_back(clusterIndex_t1[i].size());
        }
        
    /*Analize*/
        //compute the intercluster flow by intersection
        std::vector<std::vector<int>> intersectionList;//the list of minimal intersection element;every vector in list is an intersection and every inner element is index of monomer
        std::vector<std::vector<int>> lossSectionIndexList(nclust0);//lost sections at t0
        std::vector<std::vector<int>> gainSectionIndexList(nclust1);//gain sections at t1
        std::vector<int> EndLossCounter(max_sections,0);
        std::vector<int> SourceGainCounter(max_sections,0);

        for (int i = 0; i < nclust0; i++)
        {
            for (int j = 0; j < nclust1; j++)
            {
                std::vector<int> tempIntersec;
                std::set_intersection(clusterIndex_t0[i].begin(), clusterIndex_t0[i].end(), clusterIndex_t1[j].begin(), clusterIndex_t1[j].end(), std::back_inserter(tempIntersec));
                if (tempIntersec.size()!=0)
                {
                    int current_index = intersectionList.size();
                    clusterIndex_t0[i].erase(std::remove_if(clusterIndex_t0[i].begin(), clusterIndex_t0[i].end(), [&tempIntersec](const int& val){ return std::find(tempIntersec.begin(), tempIntersec.end(), val) != tempIntersec.end(); }), clusterIndex_t0[i].end());
                    clusterIndex_t1[j].erase(std::remove_if(clusterIndex_t1[j].begin(), clusterIndex_t1[j].end(), [&tempIntersec](const int& val){ return std::find(tempIntersec.begin(), tempIntersec.end(), val) != tempIntersec.end(); }), clusterIndex_t1[j].end());
                    intersectionList.push_back(tempIntersec);
                    lossSectionIndexList[i].push_back(current_index);
                    gainSectionIndexList[j].push_back(current_index);
                }
            }
        }
        //compute sections lost's contribution to alpha and beta
        for (int i = 0; i < nclust0; i++)
        {
            int nSectionLoss = lossSectionIndexList[i].size();
            int nMonomerLoss = clusterIndex_t0[i].size();
            int nEndLoss = nSectionLoss + nMonomerLoss;//number of ways this cluster break into

            if (nEndLoss < max_sections)
            {
                EndLossCounter[nEndLoss] ++;
            }
            else
            {
                EndLossCounter[max_sections-1] ++;
            }

            if (nEndLoss <= 1)
            {//==0 means break to only monomers; ==1 means completely flow into a cluster
                continue;
            }

            int cSize = cSize_t0[i];
            //compute monomer loss firstly
            for (int j = 0; j < nMonomerLoss; j++)
            {
                if (nEndLoss <= 1)
                {
                    break;
                }
                //alpha[cSize][1] += 1 / clusterDensity[cSize_t0[i]];
                alpha[cSize][1] += 1;
                cSize -= 1;
                nEndLoss -= 1;
            }
            
            for (int j = 0; j < nSectionLoss; j++)//avoid the last part of aggregation?
            {
                if (nEndLoss <= 1)
                {
                    break;
                }
                
                int sectionIndex = lossSectionIndexList[i][j];
                int sectionSize = intersectionList[sectionIndex].size();
                //alpha[cSize][sectionSize] += 1 / clusterDensity[cSize_t0[i]];
                alpha[cSize][sectionSize] += 1;
                cSize -= sectionSize;

                nEndLoss -= 1;
            }
        }
        

        for (int i = 0; i < nclust1; i++)
        {
            int nSectionGain = gainSectionIndexList[i].size();
            int nMonomerGain = clusterIndex_t1[i].size();
            int nSourceGain = nSectionGain + nMonomerGain;//number of ways this cluster get from
            
            if (nSourceGain < max_sections)
            {
                SourceGainCounter[nSourceGain] ++;
            }
            else
            {
                SourceGainCounter[max_sections - 1] ++;
            }
            

            if (nSourceGain <= 1)
            {//==0 means formed from only monomers; ==1 means completely form from a section
                continue;
            }
            
            int cSize = 0;//form cluster from  0 
            double dc=0,ds=0;
            for (int j = 0; j < nSectionGain; j++)
            {
                int sectionIndex = gainSectionIndexList[i][j];
                int sectionSize = intersectionList[sectionIndex].size();
                int sectionOriginalIndex = intersectionList[sectionIndex][0];
                int sectionOriginalSize = iSize_t0[sectionOriginalIndex];
                if (j == 0)
                {
                    dc = clusterDensity[sectionOriginalSize];
                }
                else
                {
                    ds = clusterDensity[sectionOriginalSize];
                    //beta[cSize][sectionSize] += 1/(dc*ds);
                    beta[cSize][sectionSize] += 1;
                    dc = ds;
                }
                cSize += sectionSize;
            }
            //analize monomers' contribution to beta, after the set of sections
            for (int j = 0; j < nMonomerGain; j++)
            {
                if (cSize == 0)
                {//no section flow into this cluster
                    dc = clusterDensity[1];
                }
                else
                {
                    //beta[cSize][1] += 1 / (dc*clusterDensity[1]);
                    beta[cSize][1] += 1;
                }
                
                cSize += 1;
            }
        }
        
        
        std::string alphaFilename = "flow/conf" + std::to_string(fileNameSerial) + "_alpha.dat";
        std::string betaFilename = "flow/conf" + std::to_string(fileNameSerial) + "_beta.dat";
        std::string lossFilename = "log/conf" + std::to_string(fileNameSerial) + "_loss.dat";
        std::string gainFilename = "log/conf" + std::to_string(fileNameSerial) + "_gain.dat";
        
        std::ofstream alpha_file(alphaFilename);
        
        std::ofstream beta_file(betaFilename );
        for (int i = 1; i < beta.size(); i++)
        {
            for (int j = 1; j < beta[i].size(); j++)
            {

                beta_file << beta[i][j];
                
                alpha_file << alpha[i][j];
                
                if (j != beta[i].size() - 1)
                {
                    beta_file << ",";
                    alpha_file << ",";
                }
            }
            beta_file << "\n";
            alpha_file << "\n";
        }
        beta_file.close();
        alpha_file.close();

        std::ofstream loss_file(lossFilename);

        std::ofstream gain_file(gainFilename);
        for (int i = 0; i < max_sections; i++)
        {
            loss_file<<EndLossCounter[i]<<"\n";
            gain_file<<SourceGainCounter[i]<<"\n";
        }
        loss_file.close();
        gain_file.close();

        std::ofstream logCt0("nclust0.log");
        std::ofstream logCt1("nclust1.log");
        logCt0 << nclust0 <<std::endl;
        logCt1 << nclust1 <<std::endl;
        logCt0.close();
        logCt1.close();
        

    }
    finish = clock(); 	
    printf( "Tracing cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
}

void outputAvgSize(int init_ser,int end_ser,int n_all,int maxTypeNum)
{
    clock_t start, finish;
    start = clock(); 

    double inverseV = 1.5625e-5;
    double tau = 0.2;
    double eps = 1e-12;

    vector<double> avgSizeList;

    for (size_t fileNameSerial = init_ser; fileNameSerial < end_ser; fileNameSerial++)
    {
        std::vector<std::vector<int>> clusterIndex_t0;
        std::vector<int> iSize_t0(n_all, 1);
    /*Read File*/
        std::string clustIDFilename0 = "Index/conf" + std::to_string(fileNameSerial) + "_id.dat";
        std::ifstream file;
        std::string str;
        file.open(clustIDFilename0);
        while (std::getline(file, str)) {
            std::vector<int> row;
            std::stringstream ss(str);
            int value;
            while (ss >> value) {
                row.push_back(value);
                if (ss.peek() == ',') ss.ignore();
            }
            clusterIndex_t0.push_back(row);
        }
        file.close();
    /*Count avg size*/
        int nclust0 = clusterIndex_t0.size();
        int totalElementsNum = 0;
        for (int i = 0; i < nclust0; i++)
        {
            int current_size = clusterIndex_t0[i].size();
            totalElementsNum += current_size;
        }
        double currAvgSize = (1.00*totalElementsNum) / (1.00*nclust0);
        avgSizeList.push_back(currAvgSize);
    }

    std::ofstream output;
    output.open("avgSize.log");
    for (int i = 0; i < avgSizeList.size(); i++)
    {
        output<<avgSizeList[i]<<"\n";
    }
    

    finish = clock(); 	
    printf( "statics cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
}

void inhomoAnaly(int edgeNumber)
{
    int a = edgeNumber;
    int Natom = 32;
    int Nmole = 400;
    int nline = Natom * Nmole;
    int init_fileNameIndex = 0;
    int end_fileNameIndex = 5000;
    int deltaT = 100;

    for (int i_filename = init_fileNameIndex; i_filename < end_fileNameIndex; i_filename += deltaT)
    {
        /*Read file*/
        std::string filename = "conf/conf" + std::to_string(i_filename) + ".gro";
        std::string distLogFileName = "iha_log/iha" + std::to_string(i_filename) + ".dat";
        std::vector<double> borderXYZ;
        std::vector<Point> points;
        readFile(filename,points,nline,borderXYZ);
        double L = borderXYZ[0];

        // 初始化每个小方块的点计数
        std::vector<std::vector<std::vector<int>>> count(a, std::vector<std::vector<int>>(a, std::vector<int>(a, 0)));

        // 统计每个小方块含有的点个数
        double step = L / a;
        for (const auto& p : points) {
            int ix = std::min(static_cast<int>(p.x / step), a - 1);
            int iy = std::min(static_cast<int>(p.y / step), a - 1);
            int iz = std::min(static_cast<int>(p.z / step), a - 1);
            count[ix][iy][iz]++;
        }

        // 将count[i][j][k]的数值从大到小排序，保存在一个一维数组中
        std::vector<int> count_values;
        for (int i = 0; i < a; ++i) {
            for (int j = 0; j < a; ++j) {
                for (int k = 0; k < a; ++k) {
                    count_values.push_back(count[i][j][k]);
                }
            }
        }
        std::sort(count_values.begin(), count_values.end(), std::greater<int>());

        // 将排序后的数组以每行一个值的形式写入iha_log/iha_sort1.txt
        std::ofstream sortfile(distLogFileName);
        for (int i = 0; i < count_values.size(); ++i) {
            sortfile << count_values[i] << "\n";
        }
        sortfile.close();
    }
}