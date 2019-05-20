#ifndef ADJLIST_H
#define ADJLIST_H
#include<vector>
#include<string>
const int n=62;
using namespace std;
class AdjList
{
    public:
        AdjList();
        virtual ~AdjList();
        vector < pair<int ,string> > Nodes_Labels(); // returns names of nodes
        vector < pair<int ,double> >* Main_AdjList(); // Returns adjacency list
        vector< double > Threshold();   // returns threshold of each nodes
    protected:
        //int n;
    private:
};

#endif // ADJLIST_H
