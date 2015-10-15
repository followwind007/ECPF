//============================================================================
// Name        : CSize.h
// Author      : Zhongyang Zuo
// Version     :
// Copyright   :
// Description : Ansi-style
//============================================================================


#include <stdio.h>
#include <iostream>
#include <map>

using namespace std;

class FamilyFilter {
protected:

    string clusterFileLocation;
    string filteredFileLocation;

    map<string, string> familySequence;

public:

    FamilyFilter(string cFL, string fFL);
    void filterClusterFile();

};

FamilyFilter::FamilyFilter(string cFL, string fFL){
    clusterFileLocation = cFL;
    filteredFileLocation = fFL;
}
