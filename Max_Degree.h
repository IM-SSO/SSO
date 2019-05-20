#ifndef MAX_DEGREE_H
#define MAX_DEGREE_H
#include<vector>
#include "AdjList.h"
#include "Influence_IC.h"
#include "Influence_LT.h"
using namespace std;

class Max_Degree
{
    public:
        Max_Degree();
        virtual ~Max_Degree();
        vector<int>gen_seeds(int k);
        bool active[n];
        pair<int , int> calc_Influence(vector<int>& seed_set);

    protected:

    private:
};

#endif // MAX_DEGREE_H
