#include "Max_Degree.h"

Max_Degree::Max_Degree()
{
    //ctor
}

Max_Degree::~Max_Degree()
{
    //dtor
}

vector<int> Max_Degree::gen_seeds(int k){
    for(int i=0;i<n;i++){
        active[i]=false;
    }
    AdjList a;
    vector < pair<int ,double> > *adj=a.Main_AdjList();
    vector<int>seed_set;
    for(int i=0;i<k;i++){
        int maxs=0, maxj;

        for(int j=0;j<n;j++){
            if(active[j]==true)continue;
            if(maxs<adj[j].size()){
                maxs=adj[j].size();maxj=j;
            }
            }
            seed_set.push_back(maxj);
            active[maxj]=true;

    }
    return seed_set;

}
pair<int, int> Max_Degree::calc_Influence(vector<int>&seed_set){
    pair<int, int>p;
    int a,b;
    Influence_LT obLT;

    //cout<<"Calculating influence for LT model :";
   a=obLT.Influence_set(seed_set);
    //cout<<"a ::"<<a<<endl;
     p.first=a;
    Influence_IC obIC;

    //cout<<"Calculating influence for IC model :";
    b=obIC.average_Influence_IC(seed_set, 50);//cout<<"b ::"<<b<<endl;


    p.second=b;
    return p;

}
