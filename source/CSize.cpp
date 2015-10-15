//============================================================================
// Name        : CSize.cpp
// Author      : Zhongyang Zuo
// Version     :
// Copyright   :
// Description : Ansi-style
//============================================================================


#include "CSize.h"
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>

using namespace std;

void FamilyFilter::filterClusterFile(){
    ifstream clusterStream;
    clusterStream.open(clusterFileLocation);

    if (!clusterStream.good()) {
        cerr<<"Error in opening cluster file: "<<clusterFileLocation<<endl;
    }

    string lineInfo;
    string familyID, sequenceID;
    getline(clusterStream, lineInfo);

    while (lineInfo.length()) {
        istringstream lineStream(lineInfo);
        lineStream >> familyID;
        lineStream >> sequenceID;
        familySequence[sequenceID] = familyID;

        lineInfo.clear();
        getline(clusterStream, lineInfo);

    }
    cout<<"finish in reading cluster file: "<<clusterFileLocation<<endl;

    cout<<"start to write file."<<endl;

    ofstream fileteredCluterStream;
    fileteredCluterStream.open(filteredFileLocation);

    if (!fileteredCluterStream.good()) {
        cerr<<"Error in opening filtered file: "<<filteredFileLocation<<endl;
    }

    map<string, string>::iterator fC_it;

    if (familySequence.empty()) {
        cerr<<"nothing in map."<<endl;
    }

    for (fC_it = familySequence.begin(); fC_it != familySequence.end(); fC_it++) {

        if (fC_it -> first != "ACEAZ_1_PE1875") {
            fileteredCluterStream<<fC_it -> second<<"\t"<<fC_it -> first<<"\n";
        }
        else{
            break;
        }
    }
    cout<<"success in writting filtered cluster file: "<<filteredFileLocation<<endl;

}
int main(){
    FamilyFilter ff("/Users/zuozhongyang/Wkspcs/seq.fnodes",
                    "/Users/zuozhongyang/Wkspcs/filteredseq.fnodes");
    ff.filterClusterFile();
}










