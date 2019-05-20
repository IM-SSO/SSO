#ifndef INFLUENCE_LT_H
#define INFLUENCE_LT_H
#include "AdjList.h"

using namespace std;
#include<queue>
#include<vector>


class Influence_LT
{
    public:
        Influence_LT();
        virtual ~Influence_LT();
        void init_applied(double *applied);
        void init_visited(bool *visited);

        double applied[n];
        bool visited[n];
        int Influence_set(vector<int>&v);
        int Influence_elem(int root);
    protected:

    private:
};

#endif // INFLUENCE_LT_H
