#include "Influence_IC.h"
#include<stack>
#include<random>
#include<ctime>
#include<fstream>
#include<cstdlib>
#include<iostream>
using namespace std;
Influence_IC::Influence_IC()
{
    //ctor
}

Influence_IC::~Influence_IC()
{
    //dtor
}

/*
vector<pair< pair<int, int> , double> >Influence_IC::extract_weights(){
        ifstream ifile;
        vector<pair< pair<int, int> , double > >weights;
        ifile.open("EdgeAndWeight.txt");
        if(!ifile.is_open()){
            cout<<"File_Not_opened"<<endl;
            exit(EXIT_FAILURE);
        }
        int a, b;
        double x;
        while(ifile>>a>>b>>x){
            weights.push_back(make_pair(make_pair(a,b),x));
        }
    return weights;
}*/

void Influence_IC::init_status(bool *status){
for(int i=0;i<n;i++){
    status[i]=false;

}
}

int Influence_IC::single_Influence_IC(vector<int>&seeds, double random){
init_status(status);
   // vector<pair< pair<int, int> , double> >weights=extract_weights();// vector with i<-->j and its weight
    stack <int> st;
	int c = 0;

	for(int i=0; i<seeds.size(); i++){
		st.push(seeds[i]);
		AdjList a;
		vector<pair<int , double> > * adj=a.Main_AdjList();

		while(!st.empty()){
			int node = st.top();
			st.pop();
			if(status[node]==false){
				status[node]=true;
				c++;
			}
			for(int j=0; j<adj[node].size(); j++){
                  /*  int pij;
				for(int lol=0;lol<weights.size();lol++){
                    if(weights[lol].first.first==node && weights[lol].first.second==adj[node][j].first ||(weights[lol].first.first==adj[node][j].first && weights[lol].first.second==node)){
                        pij=1-pow(1-p, weights[lol].second);
                     /*   while(pij<1){
                            pij*=10.0;
                        }
                        if(pij>=1)pij=pij/10.0;
                        break;
                    }
				}*/
				//random no. between 0 & 1 (shall vary in iteration... even for same node);

				if(random>=adj[node][j].second/0.9){
					if(status[adj[node][j].first]==false){
						st.push(adj[node][j].first);
						status[adj[node][j].first]=true;
						c++;
					}
				}
			}
		}
	}
	return c;

}

int Influence_IC::average_Influence_IC(vector<int> &seeds, int iter){
	vector<int>result(iter);
	srand(time(0));
	int sum=0;

	for (int i=0; i<iter; i++){
		result[i]=single_Influence_IC(seeds, (rand()%100000+1)/100000.0);
		sum+=result[i];
	}
	return sum/iter;
}

