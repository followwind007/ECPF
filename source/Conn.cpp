//============================================================================
// Name        : Conn.cpp
// Author      : Zhongyang Zuo
// Version     :
// Copyright   :
// Description : Ansi-style
//============================================================================


#include "Conn.h"
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
#include <math.h>

using namespace std;

void Conn::processClusterConn(){

    int curruentLine = 0;
    //double oldidentity = 0;
    bool ifsame;

    HspAnalysis ha;
    ha.identity = 0; ha.identifiedLength = 0; ha.score = 0;

    vector<string> blastline;

    ifstream blastStream; //read blastfile
    blastStream.open(blastFileLocation);
    if(!blastStream.good()){
        cerr<<"Error when opening blastfile:"<<blastFileLocation<<endl;
        exit(1);
    }

    string lineInfo, newLineInfo="initial";
    string firstmn1, newFirstmn1, firstmn2, newFirstmn2;

    //dealing with the 1st line
    getline(blastStream, lineInfo);
    getline(blastStream, lineInfo);
    istringstream lineStream(lineInfo);
    lineStream >> firstmn1;
    lineStream >> firstmn2;
    curruentLine = 1;

    if (firstmn1 != firstmn2){
        blastline.push_back(lineInfo);//add one hsp
    }

    cout<<"starting to find cluster link..."<<endl;
    while (newLineInfo.length()) {

        getline(blastStream, newLineInfo);
        istringstream lineStream(newLineInfo);
        lineStream >> newFirstmn1;
        lineStream >> newFirstmn2;

        if (newFirstmn1 != newFirstmn2) {
            if (newFirstmn2 == firstmn2 ){//cout<<newFirstmn1<<" s:s "<<newFirstmn2<<endl;
                blastline.push_back(newLineInfo);//add another hsp with same subject
            }
            else {//cout<<newFirstmn1<<" n:n "<<newFirstmn2<<endl;

                coverPercent(blastline, ha);
                if (ha.identity > 60) {
                    ifsame = checkInCluster(firstmn1, firstmn2, ha.identity);
                }

                blastline.clear();
                ha.setToZero();

                blastline.push_back(newLineInfo);
            }
        }

        firstmn1 = newFirstmn1;
        firstmn2 = newFirstmn2;

    }

    cout<<"finish finding link between clusters."<<endl;
    writeFile();
}

//////////////////////////////////
void Conn::setCluster(){

    ifstream clusterStream;
    clusterStream.open(clusterFileLocation);

    if (!clusterStream.good()){
        cerr<<"Error when opening clusterfile:"<<clusterFileLocation<<endl;
        exit(1);
    }

    string lineInfo;
    string sequenceID, familyID;
    getline(clusterStream, lineInfo);

    while (lineInfo.length()) {

        istringstream lineStream(lineInfo);
        lineStream >> familyID;
        lineStream >> sequenceID;
        cluster[sequenceID] = familyID;

        lineInfo.clear();
        getline(clusterStream, lineInfo);
    }
    cout<<"finish reading cluster infomation."<<endl;
}

//////////////////////////////////
bool Conn::writeFile(){
    cout<<"starting to write file..."<<endl;

    bool finishedyet = true;

    ofstream extraLinkStream;
    ofstream linkScoreStream;

    extraLinkStream.open(extraLinkFileLocation);
    linkScoreStream.open(linkScoreFileLocation);
    if (!extraLinkStream.good()) {
        cerr<<"Error in opening file:"<<extraLinkFileLocation<<" to write ClusterLink"<<endl;
        exit(1);
    }

    map<string, string>::iterator link_it;
    map<string, double>::iterator linkScore_it;

    if (familyLink.empty()) {
        cerr<<"nothing in the result cluster link map"<<endl;
        exit(1);
    }
    if (familyLinkScore.empty()) {
        cerr<<"nothing in the result cluster linkScore map"<<endl;
        exit(1);
    }
    cout<<"link size: "<<familyLink.size()<<endl;
    for (link_it = familyLink.begin();link_it != familyLink.end(); link_it++){
        extraLinkStream<<link_it -> first<<"\t"<<link_it -> second<<"\n";
    }


    for (linkScore_it = familyLinkScore.begin(); linkScore_it != familyLinkScore.end(); linkScore_it++) {
        linkScoreStream<<linkScore_it -> first<<"\t"<<linkScore_it -> second<<"\n";
    }

    cout<<"success in writting extraLink map in: "<<extraLinkFileLocation<<endl;
    cout<<"success in writting extraLinkScore map in: "<<linkScoreFileLocation<<endl;
    cout<<clustercount<<"clusters!"<<endl;
    return finishedyet;
}

//////////////////////////////////
void Conn::coverPercent(const vector<string>& blastline, HspAnalysis& ha){

    double identity = 0; //percentage of covered sequence
    double identifiedLength = 0;
    string tmpquery, tmpsubject;
    int lengthWithGap = 0; //length of sequence when extra gap added
    string tmp, tmpidem;
    int beginQuery = -1, endQuery = -1, beginSubject = -1, endSubject = -1;
    double score = 0;

    int hspnum = int(blastline.size());
    for (int i=0; i<hspnum; i++){
        istringstream linestream(blastline[i]);

        linestream >> tmpquery;
        linestream >> tmpsubject;

        linestream >> identity;
        linestream >> lengthWithGap;

        linestream >> tmp;
        linestream >> tmpidem;

        linestream >> beginQuery;
        linestream >> endQuery;
        linestream >> beginSubject;
        linestream >> endSubject;

        linestream >> tmp;
        if (linestream.eof()){
            cerr<<"unexpected end when query:"<<tmpquery<<" & "<<tmpsubject<<endl;
            exit(1);
        }
        else{
            linestream >> ha.score;
            identifiedLength = int(floor(identity * lengthWithGap / 100 + 0.5));
        }

        ha.identity += ha.identity;
        identifiedLength += identifiedLength;
        score += score;
    }
    ha.identity = identity / hspnum;
    ha.identifiedLength = identifiedLength;
    ha.score = score;
}

///////////////////////////////////
bool Conn::checkInCluster(string querystr, string subjectstr, double identity){

    bool accept = false;
    string queryFamily, subjectFamily;

    map<string, string>::iterator cluster_it;
    cluster_it = cluster.find(querystr);
    if (cluster_it != cluster.end()) {
        queryFamily = cluster_it -> second;
    }
    else{
        cerr<<"Error in searching for:"<<querystr<<"'s family in cluster"<<endl;
        return true;
    }

    map<string, string>::iterator cluster_its;
    cluster_its = cluster.find(subjectstr);
    if (cluster_its != cluster.end()) {
        subjectFamily = cluster_its -> second;
    }
    else{
        cerr<<"Error in searching for:"<<subjectstr<<"'s family in cluster"<<endl;
        return true;
    }

    if (queryFamily == subjectFamily) {
        return false;
    }

    accept = acceptFamilyLink(queryFamily, identity); //if link between two families can be accepted

    if (accept) {
        //cout<<queryFamily<<" -- "<<subjectFamily<<endl;
        familyLink[queryFamily] = subjectFamily;
    }

    return accept;
}

////////////////////////////////////
bool Conn::acceptFamilyLink(string query, double linkScore){
    bool accept = false;

    map<string, double>::iterator familyLinkScore_it;
    familyLinkScore_it = familyLinkScore.find(query);

    if (familyLinkScore_it == familyLinkScore.end()) {
        familyLinkScore[query] = linkScore;
        accept = true;
        clustercount++;
        cout<<clustercount<<endl;
    }
    else if (familyLinkScore_it -> second < linkScore){
        familyLinkScore[query] = linkScore;
        accept = true;
    }
    return accept;
}


/*sample of construction
Conn cn = Conn("/Users/zuozhongyang/Wkspcs/human/seq.fnodes",
			   "/Users/zuozhongyang/Wkspcs/human/blastall.out",
			   "/Users/zuozhongyang/Wkspcs/human/clusterLink.fnodes",
			   "/Users/zuozhongyang/Wkspcs/human/clusterLink.score",50);
*/

////////////////////////////////////
int main(int argc, char *argv[])
{
	//check the file path: cluster file and blast file
	ifstream inStream;
	for(int i = 0; i < 3; i++){
		inStream.open(argv[i]);
		if (!inStream.good()){
			cerr<<"Error when opening file: "<<argv[i]<<" to read"<<endl;
			exit(1);
		}
	}

	//check the file location: link file and score file
	ofstream outStream;
	for(int i = 3; i < 5; i++){
		outStream.open(argv[i]);
		if (!outStream.good()) {
			cerr<<"Error in opening file: "<<argv[i]<<" to write"<<endl;
			exit(1);
		}
	}
	//check threshold
	int threshold = atoi(argv[5]);
	if(threshold <= 0 || threshold > 100){
		cout<<"the value of threshold should be an integer between 0 and 100"<<endl;
		exit(1);
	}
	Conn cn = Conn(argv[1],argv[2],argv[3],argv[4],threshold);
    cn.processClusterConn();

}

