
#ifndef __test__Conn__
#define __test__Conn__

#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>
using namespace std;

#endif /* defined(__test__Conn__) */

struct HspAnalysis{
    double identity;
    double score;
    int identifiedLength;
    
    
    inline void setToZero(){
        identity = 0;
        identifiedLength = 0;
        score = 0;
    }
};

class Conn {
protected:
    string clusterFileLocation;
    string blastFileLocation;
    string extraLinkFileLocation;
    string linkScoreFileLocation;
    int threshold;
    int clustercount;
    
    map<string, string> cluster; //cluster infomation
    //map<string, string> extraLink; //to save the result
    map<string, string> familyLink; //link infomation between families
    map<string, double> familyLinkScore; //score that means the connection quality between two families
    
public:
    Conn(string cFL, string bFL, string eLFL, string lSFL, int th);
    
    void processClusterConn();
    
    void setCluster();
    
    bool writeFile();
    
    void coverPercent(const vector<string>& blastline, HspAnalysis& ha);
    
    bool checkInCluster(string querystr, string subjectstr, double identity);
    
    bool acceptFamilyLink(string query, double linkScore);
    
    virtual ~Conn(){};
};

Conn::Conn(string cFL, string bFL, string eLFL, string lSFL, int th){
    clusterFileLocation = cFL;
    blastFileLocation = bFL;
    extraLinkFileLocation = eLFL;
    linkScoreFileLocation = lSFL;
    threshold = th;
    clustercount = 0;
    setCluster();
}

