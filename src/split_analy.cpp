#include "split_analy.h"

void intersecTraceBack()
{
    clock_t start, finish;
    start = clock(); 
    int init_ser = 0;
    int end_ser = 200;
    //int n_all = 800;
    int n_all = 51200;
    //int maxTypeNum = 300;
    int maxTypeNum = 8800;
    double inverseV = 1.5625e-5;
    double tau = 0.2;
    double eps = 1e-9;
    int delta_t = 1;
    for (size_t fileNameSerial = init_ser; fileNameSerial < end_ser; fileNameSerial += delta_t)
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
        std::string clustIDFilename1 = "Index/conf" + std::to_string(fileNameSerial+delta_t) + "_id.dat";
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
        std::vector<int> probabilityFilter(maxTypeNum,0);
        int max_sections = 1000;
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
        /*compute section density*/
        std::vector<double> sectionDensity(maxTypeNum,0.0);
        int nMonomerS = n_all;
        for (int i = 0; i < intersectionList.size(); i++)
        {
            int current_size = intersectionList[i].size();
            if (current_size < maxTypeNum)
            {
                sectionDensity[current_size] += inverseV;
            }
            nMonomerS -= current_size;
        }
        sectionDensity[1] += nMonomerS * inverseV;

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
            double dc = 0.0;
            //compute monomer loss firstly
            for (int j = 0; j < nMonomerLoss; j++)
            {
                if (nEndLoss <= 1)
                {
                    break;
                }
                
                //alpha[cSize][1] += 1 / dc;
                alpha[cSize][1] += 1 / clusterDensity[cSize_t0[i]];
                //alpha[cSize][1] += 1;
                probabilityFilter[cSize] ++;
                cSize -= 1;
                nEndLoss -= 1;
            }
            
            for (int j = 0; j < nSectionLoss; j++)//avoid the last part of aggregation?
            {
                if (nEndLoss <= 1)
                {
                    break;
                }
                /*
                if (cSize == cSize_t0[i])
                {
                    dc = clusterDensity[cSize_t0[i]];
                }
                else if (sectionDensity[cSize] < eps)
                {
                    dc = clusterDensity[cSize_t0[i]];
                }
                else
                {
                    dc = sectionDensity[cSize];
                }*/
                
                int sectionIndex = lossSectionIndexList[i][j];
                int sectionSize = intersectionList[sectionIndex].size();
                alpha[cSize][sectionSize] += 1 / clusterDensity[cSize_t0[i]];//v3
                //alpha[cSize_t0[i]][sectionSize] += 1 / clusterDensity[cSize_t0[i]];//v5
                //alpha[cSize][sectionSize] += 1 / dc;v4.5
                //alpha[cSize][sectionSize] += 1;
                probabilityFilter[cSize] ++;
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
                    //dc = clusterDensity[sectionOriginalSize];
                    dc = sectionDensity[sectionSize];
                }
                else
                {
                    //ds = clusterDensity[sectionOriginalSize];
                    ds = sectionDensity[sectionSize];
                    beta[cSize][sectionSize] += 1/(dc*ds);
                    //beta[cSize][sectionSize] += 1;
                    dc = ds;
                    probabilityFilter[cSize] ++;
                }
                cSize += sectionSize;
            }
            //analize monomers' contribution to beta, after the set of sections
            for (int j = 0; j < nMonomerGain; j++)
            {
                if (cSize == 0)
                {//no section flow into this cluster
                    //dc = clusterDensity[1];
                    dc = sectionDensity[1];
                }
                else
                {
                    //beta[cSize][1] += 1 / (dc*clusterDensity[1]);
                    beta[cSize][1] += 1 / (dc*sectionDensity[1]);
                    //beta[cSize][1] += 1;
                    probabilityFilter[cSize] ++;
                }
                
                cSize += 1;
            }
        }
        
        
        std::string alphaFilename = "rate/conf" + std::to_string(fileNameSerial) + "_alpha.dat";
        std::string betaFilename = "rate/conf" + std::to_string(fileNameSerial) + "_beta.dat";
        //std::string alphaFilename = "rate_800/conf" + std::to_string(fileNameSerial) + "_alpha.dat";
        //std::string betaFilename = "rate_800/conf" + std::to_string(fileNameSerial) + "_beta.dat";
        std::string lossFilename = "log/conf" + std::to_string(fileNameSerial) + "_loss.dat";
        std::string gainFilename = "log/conf" + std::to_string(fileNameSerial) + "_gain.dat";
        std::ofstream alpha_file(alphaFilename);
        
        std::ofstream beta_file(betaFilename );
        for (int i = 1; i < beta.size(); i++)
        {
            for (int j = 1; j < beta[i].size(); j++)
            {
                if (probabilityFilter[i]*probabilityFilter[j] == 0)
                {
                    //beta_file<<"x";
                }
                else
                {
                    //beta_file << beta[i][j];
                    //beta_file << inverseV*beta[i][j]/tau;
                    beta_file<<i<<" "<<j<<" "<<inverseV*beta[i][j]/tau<<"\n";
                }

                if(probabilityFilter[i] == 0)
                {
                    //alpha_file<<"x";
                }
                else if(i > j)
                {
                    //alpha_file << alpha[i][j];
                    //alpha_file << inverseV*alpha[i][j]/tau;
                    alpha_file<<i<<" "<<j<<" "<< inverseV*alpha[i][j]/tau<<"\n";
                }

                if (j != beta[i].size() - 1)
                {
                    //beta_file << ",";
                    //alpha_file << ",";
                }
            }
            //beta_file << "\n";
            //alpha_file << "\n";
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
    }
    finish = clock(); 	
    printf( "Tracing cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
}

void v7_intersecTraceBack(int init_ser,int end_ser)
{
    clock_t start, finish;
    start = clock(); 
    //int n_all = 800;
    int n_all = 400;
    //int maxTypeNum = 40;
    int maxTypeNum = 400;
    double inverseV = 6.714e-4;
    double tau = 0.02;
    double eps = 1e-9;
    int delta_t = 1;
    for (size_t fileNameSerial = init_ser; fileNameSerial < end_ser; fileNameSerial += delta_t)
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
        std::string clustIDFilename1 = "Index/conf" + std::to_string(fileNameSerial+delta_t) + "_id.dat";
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
        std::vector<int> probabilityFilter(maxTypeNum,0);
        int max_sections = 1000;
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
        /*compute section density*/
        std::vector<double> sectionDensity(maxTypeNum,0.0);
        int nMonomerS = n_all;
        for (int i = 0; i < intersectionList.size(); i++)
        {
            int current_size = intersectionList[i].size();
            if (current_size < maxTypeNum)
            {
                sectionDensity[current_size] += inverseV;
            }
            nMonomerS -= current_size;
        }
        sectionDensity[1] += nMonomerS * inverseV;

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
            double inv_dc = 1 / clusterDensity[cSize_t0[i]];
            //compute monomer loss firstly
            for (int j = 0; j < nMonomerLoss; j++)
            {
                alpha[cSize][1] += inv_dc * (1.0/double(cSize));//v7.6
                probabilityFilter[cSize] ++;

                //nEndLoss -= 1;
            }
            
            for (int j = 0; j < nSectionLoss; j++)//avoid the last part of aggregation?
            {
                int sectionIndex = lossSectionIndexList[i][j];
                int sectionSize = intersectionList[sectionIndex].size();
                alpha[cSize][sectionSize] += inv_dc * (double(sectionSize)/double(cSize));//v6//v7.6

                probabilityFilter[cSize] ++;

                //nEndLoss -= 1;
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
            //analize monomers' contribution to beta, before the set of sections
            for (int j = 0; j < nMonomerGain; j++)
            {
                if (cSize == 0)
                {//no section flow into this cluster
                    //dc = clusterDensity[1];
                    dc = sectionDensity[1];
                }
                else
                {
                    //beta[cSize][1] += 1 / (dc*clusterDensity[1]);
                    beta[cSize][1] += 1 / (dc*sectionDensity[1]);
                    //beta[cSize][1] += 1;
                    probabilityFilter[cSize] ++;
                }
                
                cSize += 1;
            }
            for (int j = 0; j < nSectionGain; j++)
            {
                int sectionIndex = gainSectionIndexList[i][j];
                int sectionSize = intersectionList[sectionIndex].size();
                int sectionOriginalIndex = intersectionList[sectionIndex][0];
                int sectionOriginalSize = iSize_t0[sectionOriginalIndex];
                if (j == 0)
                {
                    //dc = clusterDensity[sectionOriginalSize];
                    dc = sectionDensity[sectionSize];
                }
                else
                {
                    //ds = clusterDensity[sectionOriginalSize];
                    ds = sectionDensity[sectionSize];
                    beta[cSize][sectionSize] += 1/(dc*ds);
                    //beta[cSize][sectionSize] += 1;
                    dc = ds;
                    probabilityFilter[cSize] ++;
                    probabilityFilter[sectionSize] ++;//added in v7.5 at 2024/7/13 20:12
                }
                cSize += sectionSize;
            }
            
        }
        std::string alphaFilename = "rate/conf" + std::to_string(fileNameSerial) + "_alpha.dat";
        std::string betaFilename = "rate/conf" + std::to_string(fileNameSerial) + "_beta.dat";
        //std::string alphaFilename = "rate_800/conf" + std::to_string(fileNameSerial) + "_alpha.dat";
        //std::string betaFilename = "rate_800/conf" + std::to_string(fileNameSerial) + "_beta.dat";
        std::string lossFilename = "log/conf" + std::to_string(fileNameSerial) + "_loss.dat";
        std::string gainFilename = "log/conf" + std::to_string(fileNameSerial) + "_gain.dat";
        std::ofstream alpha_file(alphaFilename);
        
        std::ofstream beta_file(betaFilename );
        for (int i = 1; i < beta.size(); i++)
        {
            for (int j = 1; j < beta[i].size(); j++)
            {
                if (probabilityFilter[i]*probabilityFilter[j] == 0)
                {
                    //beta_file<<"x";
                }
                else
                {
                    //beta_file << beta[i][j];
                    //beta_file << inverseV*beta[i][j]/tau;
                    beta_file<<i<<" "<<j<<" "<<inverseV*beta[i][j]/tau<<"\n";
                }

                if(probabilityFilter[i] == 0)
                {
                    //alpha_file<<"x";
                }
                else if(i > j)
                {
                    //alpha_file << alpha[i][j];
                    //alpha_file << inverseV*alpha[i][j]/tau;
                    alpha_file<<i<<" "<<j<<" "<< inverseV*alpha[i][j]/tau<<"\n";
                }

                if (j != beta[i].size() - 1)
                {
                    //beta_file << ",";
                    //alpha_file << ",";
                }
            }
            //beta_file << "\n";
            //alpha_file << "\n";
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
    }
    finish = clock(); 	
    printf( "Tracing cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
}

void v6_intersecTraceBack(int init_ser,int end_ser)
{
    clock_t start, finish;
    start = clock(); 
    int n_all = 400;
    //int maxTypeNum = 40;
    int maxTypeNum = 400;
    double inverseV = 6.714e-4;
    double tau = 0.02;
    double eps = 1e-9;
    int delta_t = 1;
    for (size_t fileNameSerial = init_ser; fileNameSerial < end_ser; fileNameSerial += delta_t)
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
        std::string clustIDFilename1 = "Index/conf" + std::to_string(fileNameSerial+delta_t) + "_id.dat";
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
        std::vector<int> probabilityFilter(maxTypeNum,0);
        int max_sections = 1000;
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
        /*compute section density*/
        std::vector<double> sectionDensity(maxTypeNum,0.0);
        int nMonomerS = n_all;
        for (int i = 0; i < intersectionList.size(); i++)
        {
            int current_size = intersectionList[i].size();
            if (current_size < maxTypeNum)
            {
                sectionDensity[current_size] += inverseV;
            }
            nMonomerS -= current_size;
        }
        sectionDensity[1] += nMonomerS * inverseV;

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
            double inv_dc = 1 / clusterDensity[cSize_t0[i]];
            //compute monomer loss firstly
            for (int j = 0; j < nMonomerLoss; j++)
            {
                alpha[cSize][1] += inv_dc*(1.0/double(cSize));
                probabilityFilter[cSize] ++;

                //nEndLoss -= 1;
            }
            
            for (int j = 0; j < nSectionLoss; j++)//avoid the last part of aggregation?
            {
                int sectionIndex = lossSectionIndexList[i][j];
                int sectionSize = intersectionList[sectionIndex].size();
                alpha[cSize][sectionSize] += inv_dc * (double(sectionSize)/double(cSize));//v6

                probabilityFilter[cSize] ++;

                //nEndLoss -= 1;
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
                    //dc = clusterDensity[sectionOriginalSize];
                    dc = sectionDensity[sectionSize];
                }
                else
                {
                    //ds = clusterDensity[sectionOriginalSize];
                    ds = sectionDensity[sectionSize];
                    beta[cSize][sectionSize] += 1/(dc*ds);
                    //beta[cSize][sectionSize] += 1;
                    dc = ds;
                    probabilityFilter[cSize] ++;
                    probabilityFilter[sectionSize] ++;
                }
                cSize += sectionSize;
            }
            //analize monomers' contribution to beta, after the set of sections
            for (int j = 0; j < nMonomerGain; j++)
            {
                if (cSize == 0)
                {//no section flow into this cluster
                    //dc = clusterDensity[1];
                    dc = sectionDensity[1];
                }
                else
                {
                    //beta[cSize][1] += 1 / (dc*clusterDensity[1]);
                    beta[cSize][1] += 1 / (dc*sectionDensity[1]);
                    //beta[cSize][1] += 1;
                    probabilityFilter[cSize] ++;
                }
                
                cSize += 1;
            }
        }
        std::string alphaFilename = "rate/conf" + std::to_string(fileNameSerial) + "_alpha.dat";
        std::string betaFilename = "rate/conf" + std::to_string(fileNameSerial) + "_beta.dat";
        //std::string alphaFilename = "rate_800/conf" + std::to_string(fileNameSerial) + "_alpha.dat";
        //std::string betaFilename = "rate_800/conf" + std::to_string(fileNameSerial) + "_beta.dat";
        std::string lossFilename = "log/conf" + std::to_string(fileNameSerial) + "_loss.dat";
        std::string gainFilename = "log/conf" + std::to_string(fileNameSerial) + "_gain.dat";
        std::ofstream alpha_file(alphaFilename);
        
        std::ofstream beta_file(betaFilename );
        for (int i = 1; i < beta.size(); i++)
        {
            for (int j = 1; j < beta[i].size(); j++)
            {
                if (probabilityFilter[i]*probabilityFilter[j] == 0)
                {
                    //beta_file<<"x";
                }
                else
                {
                    //beta_file << beta[i][j];
                    //beta_file << inverseV*beta[i][j]/tau;
                    beta_file<<i<<" "<<j<<" "<<inverseV*beta[i][j]/tau<<"\n";
                }

                if(probabilityFilter[i] == 0)
                {
                    //alpha_file<<"x";
                }
                else if(i > j)
                {
                    //alpha_file << alpha[i][j];
                    //alpha_file << inverseV*alpha[i][j]/tau;
                    alpha_file<<i<<" "<<j<<" "<< inverseV*alpha[i][j]/tau<<"\n";
                }

                if (j != beta[i].size() - 1)
                {
                    //beta_file << ",";
                    //alpha_file << ",";
                }
            }
            //beta_file << "\n";
            //alpha_file << "\n";
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
    }
    finish = clock(); 	
    printf( "Tracing cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
}

void splitComputeAverageRate(int maxTypeNum,int beginFileSer,int endFileSer,int delta_t)
{
    splitComputeAverageAlpha(maxTypeNum,beginFileSer,endFileSer,delta_t);
    splitComputeAverageBeta(maxTypeNum,beginFileSer,endFileSer,delta_t);
}

void splitComputeAverageAlpha(int maxTypeNum,int beginFileSer,int endFileSer,int delta_t)
{
    clock_t start, finish;
    start = clock(); 
    int m = maxTypeNum;
    
    std::ofstream output;
    //compute average alpha
    std::vector<std::vector<double>> M_ave_a(m, std::vector<double>(m, 0));
    std::vector<std::vector<int>> count_a(m, std::vector<int>(m, 0));
    for (int file_num = beginFileSer; file_num < endFileSer; file_num += delta_t) {
        std::ifstream file("rate/conf" + std::to_string(file_num)+"_alpha.dat");
        std::string line;
        int row_num = 0;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            int i, j;
            double x;
            if (!(iss >> i >> j >> x)) { break; } // read i, j, x
            M_ave_a[i][j] += x;
            count_a[i][j]++;
        }
        file.close();
    }
    for (int i = 1; i < m; i++) {
        for (int j = 1; j < m; j++) {
            if (count_a[i][j] != 0) {
                M_ave_a[i][j] /= count_a[i][j];
            }
        }
    }

    //output alpha
    output.open("average_alpha.dat");
    for (size_t i = 1; i < m; i++)
    {
        for (size_t j = 1; j < m; j++)
        {
            output<<M_ave_a[i][j]<<"\n";
        }
    }
    output.close(); 
    finish = clock(); 	
    printf( "Compute average alpha cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
    
}

void splitComputeAverageBeta(int maxTypeNum,int beginFileSer,int endFileSer,int delta_t)
{
    clock_t start, finish;
    start = clock(); 
    int m = maxTypeNum;
    
    std::ofstream output;
    //compute average beta
    std::vector<std::vector<double>> M_ave(m, std::vector<double>(m, 0));
    std::vector<std::vector<int>> count(m, std::vector<int>(m, 0));
    for (int file_num = beginFileSer; file_num < endFileSer; file_num += delta_t) {
        std::ifstream file("rate/conf" + std::to_string(file_num)+"_beta.dat");
        std::string line;
        int row_num = 0;
        while (std::getline(file, line))
        {
            std::istringstream iss(line);
            int i, j;
            double x;
            if (!(iss >> i >> j >> x)) { break; } // read i, j, x
            M_ave[i][j] += x;
            count[i][j]++;
        }
        file.close();
    }
    for (int i = 1; i < m; i++) {
        for (int j = 1; j < m; j++) {
            if (count[i][j] != 0) {
                M_ave[i][j] /= count[i][j];
            }
        }
    }

    //output beta
    output.open("average_beta.dat");
    for (size_t i = 1; i < m; i++)
    {
        for (size_t j = 1; j < m; j++)
        {
            output<<M_ave[i][j]<<"\n";
        }
    }
    output.close(); 
    
    finish = clock(); 	
    printf( "Compute average beta cost:%f seconds (%f hours)\n", (double)(finish - start) / CLOCKS_PER_SEC, (double)(finish - start) / CLOCKS_PER_SEC / 3600 );
    
}

void splitCoarseGrainedAverage()
{//merge some intervals to coarse-grained bigger interval
    int n_fine_intervals = 99;
    int nType = 8799;
    int cg_gap = 100;
    int n_cg_intervals = 87;
    int n_interv = n_cg_intervals + n_fine_intervals;
    
    FILE *fp;
    //compute alpha
    std::vector<std::vector<double>> alpha_origin(nType, std::vector<double>(nType, 0));
    std::vector<std::vector<double>> alpha_cg(n_interv, std::vector<double>(n_interv, 0));
    if(fp = fopen("average_alpha.dat", "r"))
	{
		for (int i = 0; i < nType; i++)
		{
			for (int j = 0; j < nType; j++)
			{
				fscanf(fp, "%lf", &alpha_origin[i][j]);
			}
		}
		fclose(fp);
	}
    std::cout<<"reading files complete"<<"\n";
    for (int i = 0; i < n_fine_intervals; i++)
    {
        for (int j = 0; j < n_fine_intervals; j++)
        {
            alpha_cg[i][j] = alpha_origin[i][j];
        }
        for (int j = n_fine_intervals; j < n_interv; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {   //sum only at j direction
                //int current_index = (j - n_fine_intervals + 1) * cg_gap + k;
                int current_index = n_fine_intervals + (j - n_fine_intervals) * cg_gap + k;
                alpha_cg[i][j] += alpha_origin[i][current_index];
            }
            alpha_cg[i][j] /= cg_gap;
        }
    }
    std::cout<<"left up complete"<<"\n";
    for (int i = n_fine_intervals; i < n_interv; i++)
    {
        for (int j = 0; j < n_fine_intervals; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {//sum only at i direction
                //int current_index = (i - n_fine_intervals + 1) * cg_gap + k;
                int current_index = n_fine_intervals + (i - n_fine_intervals) * cg_gap + k;
                alpha_cg[i][j] += alpha_origin[current_index][j];
            }
            alpha_cg[i][j] /= cg_gap;
        }

        for (int j = n_fine_intervals; j < n_interv; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {//sum at both i and j direction
                //int current_index_i = (i - n_fine_intervals + 1) * cg_gap + k;
                int current_index_i = n_fine_intervals + (i - n_fine_intervals) * cg_gap + k;
                for (int l = 0; l < cg_gap; l++)
                {
                    //int current_index_j = (j - n_fine_intervals + 1) * cg_gap + l;
                    int current_index_j = n_fine_intervals + (j - n_fine_intervals) * cg_gap + l;
                    alpha_cg[i][j] += alpha_origin[current_index_i][current_index_j];
                }
            }
            alpha_cg[i][j] /= (cg_gap*cg_gap);
        }
    }
    std::cout<<"alpha complete"<<"\n";
    //compute beta
    std::vector<std::vector<double>> beta_origin(nType, std::vector<double>(nType, 0));
    std::vector<std::vector<double>> beta_cg(n_interv, std::vector<double>(n_interv, 0));
    if(fp = fopen("average_beta.dat", "r"))
	{
		for (int i = 0; i < nType; i++)
		{
			for (int j = 0; j < nType; j++)
			{
				fscanf(fp, "%lf", &beta_origin[i][j]);
			}
		}
		fclose(fp);
	}
    for (int i = 0; i < n_fine_intervals; i++)
    {
        for (int j = 0; j < n_fine_intervals; j++)
        {
            beta_cg[i][j] = beta_origin[i][j];
        }
        for (int j = n_fine_intervals; j < n_interv; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {   //sum only at j direction
                //int current_index = (j - n_fine_intervals + 1) * cg_gap + k;
                int current_index = n_fine_intervals + (j - n_fine_intervals) * cg_gap + k;
                beta_cg[i][j] += beta_origin[i][current_index];
            }
            beta_cg[i][j] /= cg_gap;
        }
    }
    for (int i = n_fine_intervals; i < n_interv; i++)
    {
        for (int j = 0; j < n_fine_intervals; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {//sum only at i direction
                //int current_index = (i - n_fine_intervals + 1) * cg_gap + k;
                int current_index = n_fine_intervals + (i - n_fine_intervals) * cg_gap + k;
                beta_cg[i][j] += beta_origin[current_index][j];
            }
            beta_cg[i][j] /= cg_gap;
        }

        for (int j = n_fine_intervals; j < n_interv; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {//sum at both i and j direction
                //int current_index_i = (i - n_fine_intervals + 1) * cg_gap + k;
                int current_index_i = n_fine_intervals + (i - n_fine_intervals) * cg_gap + k;
                for (int l = 0; l < cg_gap; l++)
                {
                    //int current_index_j = (j - n_fine_intervals + 1) * cg_gap + l;
                    int current_index_j = n_fine_intervals + (j - n_fine_intervals) * cg_gap + l;
                    beta_cg[i][j] += beta_origin[current_index_i][current_index_j];
                }
            }
            beta_cg[i][j] /= (cg_gap*cg_gap);
        }
    }

    std::ofstream output;
    output.open("split_cg_alpha.dat");
    for (size_t i = 0; i < n_interv; i++)
    {
        for (size_t j = 0; j < n_interv; j++)
        {
            output<<alpha_cg[i][j]<<"\n";
        }
    }
    output.close();
    output.open("split_cg_beta.dat");
    for (size_t i = 0; i < n_interv; i++)
    {
        for (size_t j = 0; j < n_interv; j++)
        {
            output<<beta_cg[i][j]<<"\n";
        }
    }
    output.close(); 
}


void splitCoarseGrainedAverage_withName(std::string inFilename,std::string outFilename)
{//merge some intervals to coarse-grained bigger interval
    int n_fine_intervals = 49;
    int nType = 4399;
    int cg_gap = 50;
    int n_cg_intervals = 87;
    int n_interv = n_cg_intervals + n_fine_intervals;
    
    //FILE *fp;
    std::ifstream readin;
    std::string line;
    //compute beta
    std::vector<std::vector<double>> beta_origin(nType, std::vector<double>(nType, 0));
    std::vector<std::vector<double>> beta_cg(n_interv, std::vector<double>(n_interv, 0));

    readin.open(inFilename);
    for (int i = 0; i < nType; i++)
    {
        for (int j = 0; j < nType; j++)
        {
            std::getline(readin,line);
            beta_origin[i][j] = std::stod(line);
        }
    }
    readin.close();
    for (int i = 0; i < n_fine_intervals; i++)
    {
        for (int j = 0; j < n_fine_intervals; j++)
        {
            beta_cg[i][j] = beta_origin[i][j];
        }
        for (int j = n_fine_intervals; j < n_interv; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {   //sum only at j direction
                //int current_index = (j - n_fine_intervals + 1) * cg_gap + k;
                int current_index = n_fine_intervals + (j - n_fine_intervals) * cg_gap + k;
                beta_cg[i][j] += beta_origin[i][current_index];
            }
            beta_cg[i][j] /= cg_gap;
        }
    }
    for (int i = n_fine_intervals; i < n_interv; i++)
    {
        for (int j = 0; j < n_fine_intervals; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {//sum only at i direction
                //int current_index = (i - n_fine_intervals + 1) * cg_gap + k;
                int current_index = n_fine_intervals + (i - n_fine_intervals) * cg_gap + k;
                beta_cg[i][j] += beta_origin[current_index][j];
            }
            beta_cg[i][j] /= cg_gap;
        }

        for (int j = n_fine_intervals; j < n_interv; j++)
        {
            for (int k = 0; k < cg_gap; k++)
            {//sum at both i and j direction
                //int current_index_i = (i - n_fine_intervals + 1) * cg_gap + k;
                int current_index_i = n_fine_intervals + (i - n_fine_intervals) * cg_gap + k;
                for (int l = 0; l < cg_gap; l++)
                {
                    //int current_index_j = (j - n_fine_intervals + 1) * cg_gap + l;
                    int current_index_j = n_fine_intervals + (j - n_fine_intervals) * cg_gap + l;
                    beta_cg[i][j] += beta_origin[current_index_i][current_index_j];
                }
            }
            beta_cg[i][j] /= (cg_gap*cg_gap);
        }
    }

    std::ofstream output;
    output.open(outFilename);
    for (size_t i = 0; i < n_interv; i++)
    {
        for (size_t j = 0; j < n_interv; j++)
        {
            output<<beta_cg[i][j]<<"\n";
        }
    }
    output.close(); 
}