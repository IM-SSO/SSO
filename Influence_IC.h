#ifndef INFLUENCE_IC_H
#define INFLUENCE_IC_H
#include "AdjList.h"

class Influence_IC
{
    public:
        Influence_IC();
        virtual ~Influence_IC();
        bool status[n];
        void init_status(bool *status);
        int average_Influence_IC(vector<int>&seeds, int iter);
        int single_Influence_IC(vector<int>&seeds, double random);
        //vector<pair<pair<int, int>, double> >extract_weights();

    protected:

    private:
};

#endif // INFLUENCE_IC_H
