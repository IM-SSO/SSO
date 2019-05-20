#ifndef RANDOMALGO_H
#define RANDOMALGO_H
#include "AdjList.h"
#include "Influence_IC.h"
#include "Influence_LT.h"

class RandomAlgo
{
    public:
        RandomAlgo();
        virtual ~RandomAlgo();

        bool active[n];     //stores current status of nodes
        void init_active(bool * active);      // initializes all nodes to false
        vector<int> Gen_seeds(int k);   // Generates seed nodes
        pair<int , int> calc_Influence(vector<int>& seed_set);

    protected:

    private:
};

#endif // RANDOMALGO_H
