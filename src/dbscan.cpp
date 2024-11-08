#include "dbscan.h"

int DBSCAN::run()
{
    int clusterID = 0;
    vector<Point>::iterator iter;
    for(iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( iter->clusterID == UNCLASSIFIED )
        {
            if ( expandCluster(*iter, clusterID) != FAILURE )
            {
                clusterID += 1;
            }
        }
    }

    return clusterID;
}

int DBSCAN::expandCluster(Point point, int clusterID)
{    
    vector<int> clusterSeeds = calculateCluster(point);

    if ( clusterSeeds.size() < m_minPoints )
    {
        point.clusterID = NOISE;
        return FAILURE;
    }
    else
    {
        int index = 0, indexCorePoint = 0;
        vector<int>::iterator iterSeeds;
        for( iterSeeds = clusterSeeds.begin(); iterSeeds != clusterSeeds.end(); ++iterSeeds)
        {
            m_points.at(*iterSeeds).clusterID = clusterID;
            if (m_points.at(*iterSeeds).x == point.x && m_points.at(*iterSeeds).y == point.y && m_points.at(*iterSeeds).z == point.z )
            {
                indexCorePoint = index;
            }
            ++index;
        }
        clusterSeeds.erase(clusterSeeds.begin()+indexCorePoint);

        for( vector<int>::size_type i = 0, n = clusterSeeds.size(); i < n; ++i )
        {
            vector<int> clusterNeighors = calculateCluster(m_points.at(clusterSeeds[i]));

            if ( clusterNeighors.size() >= m_minPoints )
            {
                vector<int>::iterator iterNeighors;
                for ( iterNeighors = clusterNeighors.begin(); iterNeighors != clusterNeighors.end(); ++iterNeighors )
                {
                    if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED || m_points.at(*iterNeighors).clusterID == NOISE )
                    {
                        if ( m_points.at(*iterNeighors).clusterID == UNCLASSIFIED )
                        {
                            clusterSeeds.push_back(*iterNeighors);
                            n = clusterSeeds.size();
                        }
                        m_points.at(*iterNeighors).clusterID = clusterID;
                    }
                }
            }
        }

        return SUCCESS;
    }
}

vector<int> DBSCAN::calculateCluster(Point point)
{
    int index = 0;
    vector<Point>::iterator iter;
    vector<int> clusterIndex;
    for( iter = m_points.begin(); iter != m_points.end(); ++iter)
    {
        if ( calculateDistance(point, *iter) <= m_epsilon )
        {
            clusterIndex.push_back(index);
        }
        index++;
    }
    return clusterIndex;
}

inline double DBSCAN::calculateDistance(const Point& pointCore, const Point& pointTarget )
{
    double ceilx = m_borderList[0];
    double ceily = m_borderList[1];
    double ceilz = m_borderList[2];
    double dx,dy,dz = 0.0;
    double xc = pointCore.x;
    double xt = pointTarget.x;
    double yc = pointCore.y;
    double yt = pointTarget.y;
    double zc = pointCore.z;
    double zt = pointTarget.z;

    if (xc < m_rcut && xt >(ceilx-m_rcut))
        dx = xc - xt + ceilx;
    else if (xt < m_rcut && xc >(ceilx-m_rcut))
        dx = xt - xc + ceilx;
    else
        dx = xc - xt;
    
    if (yc < m_rcut && yt >(ceily-m_rcut))
        dy = yc - yt + ceily;
    else if (yt < m_rcut && yc >(ceily-m_rcut))
        dy = yt - yc + ceily;
    else
        dy = yc - yt;
    
    if (zc < m_rcut && zt >(ceilz-m_rcut))
        dz = zc - zt + ceilz;
    else if (zt < m_rcut && zc >(ceilz-m_rcut))
        dz = zt - zc + ceilz;
    else
        dz = zc - zt;

    return pow(dx,2)+pow(dy,2)+pow(dz,2);
    //return pow(pointCore.x - pointTarget.x,2)+pow(pointCore.y - pointTarget.y,2)+pow(pointCore.z - pointTarget.z,2);
}


