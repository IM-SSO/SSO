#ifndef GREEDY_H
#define GREEDY_H
#include "AdjList.h"
#include "Influence_IC.h"
#include "Influence_LT.h"

class Greedy
{
    public:
        Greedy();
        virtual ~Greedy();
        void init_added(bool * added);
        bool added[n];
        vector<int>seed;
        pair<vector<int>, int >final_infl(int k, int flag);
    protected:


    private:
};

#endif // GREEDY_H
