#include "Influence_LT.h"
#include<iostream   >

using namespace std;
Influence_LT::Influence_LT()
{
    //ctor
}

Influence_LT::~Influence_LT()
{
    //dtor
}



void Influence_LT::init_applied(double* applied){
            for(int i=0;i<n;i++){
                applied[i]=0;
                }
}


void Influence_LT::init_visited(bool * visited){
    for(int i=0;i<n;i++){
        visited[i]=false;;
    }
}
int Influence_LT:: Influence_elem(int root){
        int c=0;
        queue<int>q;
        q.push(root);

        if(visited[root]==false){
            visited[root]=true;
            c++;
        }
        AdjList a;

        vector<pair<int , double> >*adj=a.Main_AdjList();

        vector<double>threshold=a.Threshold();
       // cout<<"In infl_elem "<<endl;
        while(!q.empty()){
            int f=q.front();
            q.pop();
            for(int i=0;i<adj[f].size();i++){
                if(visited[adj[f][i].first]==false){
                    applied[adj[f][i].first]+=adj[f][i].second;
                    if(applied[adj[f][i].first]>threshold[adj[f][i].first -1]){
                        c++;
                        visited[adj[f][i].first]=true;
                        q.push(adj[f][i].first);
                    }
                }
            }

        }
        return c;


}
int Influence_LT::Influence_set(vector<int>&v){
    int  influence=0;

    init_applied(applied);
    init_visited(visited);

    for(int i=0;i<v.size();i++){
        int x=Influence_elem(v[i]);
        influence+=x;
        //cout<<"Element and influence :"<<v[i]<<"  "<<x<<endl;
    }

    return influence;
}


















