#include <iostream>
#include<fstream>
#include "RandomAlgo.h"
#include "Influence_LT.h"
#include "Influence_IC.h"
//const int n=62;
#include "Max_Degree.h"
#include "Greedy.h"
using namespace std;

int main()
{
    //cout<<"******************Calculating Influence through Random Algorithm for both LT and IC models******************"<<endl;

    int k;
    RandomAlgo randomobject;

    cout<<" n="<<n<<endl;

    pair<int, int>infl_res;
    vector<int>seed;
   /* ofstream rand_LT;
    ofstream rand_IC;
    rand_IC.open("random_IC.txt");
    rand_LT.open("random_LT.txt");
    for(int k=1;k<30;){
        cout<<"k = "<<k<<endl;
        seed=randomobject.Gen_seeds(k);
        for(int i=0;i<k;i++){
            cout<<seed[i]<<" ";
        }
        cout<<endl;

        infl_res.first=randomobject.calc_Influence(seed).first;
        cout<<"Through LT :"<<infl_res.first<<" ";
               infl_res.second=randomobject.calc_Influence(seed).second;

        cout<<"Through IC :"<<infl_res.second<<endl<<endl;
        rand_IC<<k<<"   "<<infl_res.second<<endl;
        rand_LT<<k<<"   "<<infl_res.first<<endl;
        k=k+2;
    }rand_IC.close();
    rand_LT.close();
     cout<<"**********************************************END******************************************************"<<endl<<endl;
    cout<<"******************Calculating Influence through Max_Degree Algorithm for both LT and IC models************************"<<endl;
    Max_Degree objMax;
    ofstream max_IC, max_LT;
    max_IC.open("maxDegree_IC.txt");
    max_LT.open("maxDegree_LT.txt");
    for(k=1;k<30;k=k+2){
        seed=objMax.gen_seeds(k);
         for(int i=0;i<k;i++){
            cout<<seed[i]<<" ";
        }
        cout<<endl;
        infl_res.first=objMax.calc_Influence(seed).first;
        cout<<"Through LT :"<<infl_res.first<<" ";
               infl_res.second=objMax.calc_Influence(seed).second;

        cout<<"Through IC :"<<infl_res.second<<endl<<endl;
        max_IC<<k<<"    "<<infl_res.second<<endl;
        max_LT<<k<<"    "<<infl_res.first<<endl;
    }
        max_IC.close();
        max_LT.close();*/

   // cout<<"**********************************************END******************************************************"<<endl<<endl;
    cout<<"******************Calculating Influence through Greedy Algorithm for both LT and IC models************************"<<endl;
    Greedy objGreedy;
    ofstream gre_IC, gre_LT;
    gre_IC.open("Greedy_IC.txt");
    gre_LT.open("Greedy_LT.txt");
    pair<vector<int>, int > p;
    for(k=1;k<30;k=k+2){
            cout<<endl;
         p=objGreedy.final_infl(k, 1);
         cout<<"k is :"<<k<<endl;
        cout<<"*****Generated seed for LT*****"<<endl;
        seed=p.first;
          cout<<"Solution set is :";
        for(int j=0;j<seed.size();j++)cout<<seed[j]<<" ";
        cout<<endl;
        cout<<"Through LT :"<<p.second<<endl;
        gre_LT<<k<<"    "<<p.second<<endl;
        cout<<"*****Generated seed for IC*****"<<endl;
        p=objGreedy.final_infl(k, 0);
        seed=p.first;
        for(int j=0;j<seed.size();j++)cout<<seed[j]<<" ";
        cout<<endl;
        cout<<"Through IC :"<<p.second<<endl;
        gre_IC<<k<<"    "<<p.second<<endl;


    }
    gre_IC.close();
    gre_LT.close();
    cout<<"**********************************************END******************************************************"<<endl<<endl;


return 0;
}
