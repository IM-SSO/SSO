#include "RandomAlgo.h"
#include<random>
#include<ctime>

#include<iostream>
RandomAlgo::RandomAlgo()
{
    //ctor
}

RandomAlgo::~RandomAlgo()
{
    //dtor
}

void RandomAlgo::init_active(bool * active){
    for(int i=0;i<n;i++){
       active[i]=false;
    }
}
vector<int>RandomAlgo::Gen_seeds(int k){
        uniform_int_distribution<int> range(0, n);
        default_random_engine generator;
        init_active(active);
        vector<int>seed_set;
        for(int i=0;i<k;i++){
                int m=range(generator);
                 if(active[m]==false){
                    seed_set.push_back(m);
                    active[m]=true;
                 }else i--;


        }
        //init_active(active);
        return seed_set;
}

pair<int, int> RandomAlgo::calc_Influence(vector<int>&seed_set){
    //int i;
    //cout<<"Influence through LT : ";
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
