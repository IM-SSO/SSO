#include "AdjList.h"
#include<fstream>
#include<vector>
#include<string>
#include<utility>
#include<cstdlib>
#include<iostream>
#include<random>
#include<ctime>

using namespace std;
AdjList::AdjList()
{
    //ctor
}

AdjList::~AdjList()
{
    //dtor
}
vector <pair<int ,string> > AdjList::Nodes_Labels(){
    vector <pair<int ,string> >Node_Names; // Vector containin names of nodes
        ifstream NodeOpen;
        NodeOpen.open("Nodes_Labels.txt");
        int num; string s;
        /*if(!NodeOpen.is_open()){
            cout<<"Can't open File_ NOdes_Labels "<<endl;
            exit(EXIT_FAILURE);
            }else cout<<"File Opened _ NOdes_Labels "<<endl;
*/
        while(NodeOpen>>num>>s){

            Node_Names.push_back(make_pair(num, s));
        }
        NodeOpen.close();
        return Node_Names;
}

vector <pair <int, double> >* AdjList:: Main_AdjList(){
        ifstream ifile;
        ifile.open("EdgeAndWeight.txt");
        /*if(!ifile.is_open()){
		cout<<"Can't open file_ EdgeAndWEight "<<endl;
		exit(EXIT_FAILURE);
	}else cout<<"File opened_EdgeAndWEight"<<endl;*/
        int n=Nodes_Labels().size();
	vector <pair <int, double> > *adj= new vector <pair <int, double> >[n];
	int x,y;double p;

	while(ifile>>x>>y>>p){

        adj[x].push_back(make_pair(y, p));
        adj[y].push_back(make_pair(x,p));

	}
	ifile.close();
    return adj;


}

vector<double> AdjList:: Threshold(){
    ifstream ifile;
    ifile.open("threshold.txt");
    double x;
  /*  if(!ifile.is_open()){
        cout<<"Can't open file_ Threshold "<<endl;
		exit(EXIT_FAILURE);
    }
    else cout<<"File Opened_Threshold "<<endl;*/
    vector<double >threshold;
    while(ifile>>x){
        threshold.push_back(x);
    }
    return threshold;
}




