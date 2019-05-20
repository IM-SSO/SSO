#include "Greedy.h"
using namespace std;
#include<iostream>
Greedy::Greedy()
{
    //ctor
}

Greedy::~Greedy()
{
    //dtor
}
void Greedy::init_added(bool *added){
    for(int i=0;i<n;i++){
        added[i]=false;
    }
}

pair<vector<int>, int > Greedy::final_infl(int k, int flag){
    vector<int>seed;
    bool added[n];
    init_added(added);
    //for(int i=0;i<n;i++)cout<<added[i]<<"le mkl"<<endl;
    int tm=0, maxi=0;
     // tm stores the maximum influence after adding each element and maxi stores its index;
    for(int i=0;i<k;i++){
//cout<<1<<endl;
        for(int j=0;j<n;j++){
            //cout<<"lol";
            if(added[j]==false){
                {
                    seed.push_back(j);
                    //cout<<"mkb";

                }
                if(flag ==1){
                        //cout<<"lol";
                    Influence_LT objLT;
                    int temp= objLT.Influence_set(seed);
                    //zzzzxczzz-cout<<temp<<"  temp"<<endl;
                    if(temp>tm){
                        tm=temp;
                        maxi=j;
                    }
                }else {
                    //cout<<"llllll"
;                    Influence_IC objIC;
                int temp=objIC.average_Influence_IC(seed, 5);
                        if(temp>tm){
                        tm=temp;
                        maxi=j;
                    }
                }
                    seed.pop_back();
            }

        }
        seed.push_back(maxi);
        added[maxi]=true;
        }
    return make_pair(seed, tm);
}
