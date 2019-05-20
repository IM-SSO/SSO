#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<cstdlib>
#include<queue>
#include<random>
#include<iomanip>
#include<climits>
#include<chrono>
#include<map>
#include<cstdio>
#include<string>
#include<cstdlib>
#include<utility>
#include<algorithm>
#define siz 430
#define ss 30
using namespace std;
vector <pair<int,double>  >adj[10001];
vector<double> threshold;
bool added[siz];
bool visited [siz];
double applied[siz];

vector<pair<double, double> > cartesian_map;
// Coordinates of nodes
vector<pair<double, double> > coordinates() {
	ifstream ifile;
	ifile.open("coordinates_final.txt");
	vector<pair<double, double> > v;
	v.resize(siz);
	if(!ifile.is_open())cout<<"Coordinates file not open "<<endl;

	int a; double b,c;
	while(ifile>>a>>b>>c){

		v[a]=make_pair(b,c);
	}
	return v;

}

void init_bool(bool * array){
	for(int i=0; i<siz; i++)
		array[i]=false;
}
void init_int(int * array){
	for(int i=0; i<siz; i++)
		array[i]=0;
}
void init_double(double *array)
{for(int i=0; i<siz; i++)
		array[i]=0;

 }
void adjlist()//passed the test
{   cout.precision(1);


ifstream ifile;
ifile.open("threshold.txt");
double p;
while(ifile>>p){
    threshold.push_back(p);
}
ifile.close();


    ifstream read;
    int nodes, edges;
    pair<int , int>x;
    read.open("EdgeAndWeight.txt");

   int a,b;
   double c;
   while(read>>a>>b>>c){
    adj[a].push_back(make_pair(b,c));
    adj[b].push_back(make_pair(a,c));

   }
    }


int  influence_elem(int root)
    {int c=0;
    queue<int> q;
    q.push(root);
	if(visited[root]==false)
	{visited[root]=true;
	 c++;
	}
    while(!q.empty())
    {
        int f=q.front();
        q.pop();
        cout<<"popped element"<<f<<endl;
        for(int i=0; i<adj[f].size(); i++)
        {if(visited[adj[f][i].first]==false)// && adj[f][i].first.second<adj[f][i].second)
            {applied[adj[f][i].first]+=adj[f][i].second;
                if(applied[adj[f][i].first]>threshold[adj[f][i].first-1]){
                      visited[adj[f][i].first]=true;
                	c++;q.push(adj[f][i].first);}

            }
        }
        }
        cout<<"influence of element"<<c<<" ";
        return c;
}



int influence_set(vector<int> &x)
{ cout<<"the set passed for check"<<endl;
for(int i=0;i<x.size();i++)
{
    cout<<x[i]<<" ";
}cout<<endl;
    int infl=0;init_bool(visited);init_double(applied);
    for(int i=0;i<x.size();i++)
    {infl+=influence_elem(x[i]);}
    cout<<"influence of set"<<infl<<endl;
    return infl;
}


    class Vibration
    {public:
        double intensity;
     };
    class Spider
    {public:
        vector<int> Position;
        int fitness;
        Vibration tar;
        vector<int> pTarPostion;
        vector<int> dimension_mask;
        vector<int> prev_move;
        int cs;
     };
    vector<int> max_degree()
{vector<pair<int,int>> indegree(430);
    for(int i=0;i<430;i++)
    {
        for(int j=0;j<adj[i].size();j++){
            indegree[adj[i][j].first].first++;
            indegree[adj[i][j].first].second=adj[i][j].first;
        }
    }
    //sort(indegree.begin(),indegree.end(),greater<pair<int,int>>());
   // cout<<"indegree sol"<<endl;
    for(int i=0;i<indegree.size();i++)
    {
        cout<<indegree[i].first<<" "<<indegree[i].second<<endl;
    }
    vector<int> solution;
    for(int i=0;i<ss;i++)
    {
        solution.push_back(indegree[i].second);
    }
    //cout<<"yahann se paaar ho gaya tha "<<endl;
return solution;

 }
double source_vibration(int influ)
{
    return log((double)1.0/(double)(influ-2)+1.0);
}
double mean(vector<double> data) {
    double sum = 0.0;
    for (int i = 0; i < data.size(); ++i) {
        sum += data[i];
    }
    return sum / data.size();
}
double std_dev(vector<double> data) {
    //cout<<".........standard deviation...........";
    double mean_val = mean(data);
    double sum = 0.0;
    for (int i = 0; i < data.size(); ++i) {
        sum += (mean_val - data[i]) * (mean_val - data[i]);
    }
    return sqrt(sum / data.size());
 }
double vibration_generation(vector<Spider> population) {
       //cout<<"...................vibration generation ................"<<endl;
    double sum = 0.0;
    vector<double> data;
    data.resize(population.size());
    for (int i = 0; i < 10; ++i) //seed set ka size
        {
        for (int j = 0; j < population.size(); ++j) {
            data[j] = (population[j].Position)[i];
        }
        sum += std_dev(data);
    }
    return sum;
      //cout<<"...................vibration generation ................"<<endl;
 }
int distance(vector<int> A,vector<int> B)
{
    double sum=0;
    for(int i=0;i<A.size();i++)
    {
	    double t;
	    t=sqrt(pow(cartesian_map[A[i]].first-cartesian_map[B[i]].first,2)+ pow(cartesian_map[A[i]].second-cartesian_map[B[i]].second,2));
	    sum+=t;
    }
    sum/=A.size();
    return (int)sum;
}
 int main()

{clock_t t1,t2;
    t1=clock();

       cartesian_map = coordinates();
    /*cout<<"hello world"<<endl;
    cout<<"enter pc & pm "<<endl;
   */ double pc, pm;
    //cin>>pc>>pm;
   pc=0.7;
   pm=0.5;
     //cout<<"::the adjacency list::"<<endl;
    adjlist();
   /* for(int i=0;i<430;i++)
    {cout<<"["<<i<<"]"<<"->";
        for(int j=0;j<adj[i].size();j++)
        {
            cout<<"["<<adj[i][j].first<<" ,"<<adj[i][j].second<<"]";
             cout<<"->";
        }//cout<<endl;
    }*/
   // cout<<"::Threshold array::"<<endl;
   /* for(int i=0;i<threshold.size();i++)
    {
        cout<<threshold[i]<< " ";
    }cout<<endl;*/
   //// cout<<"'''''calling greedy'''''"<<endl;
    vector<int> solution;
    //solution=greedy(4);
    cout<<"::Solution Seed Set::"<<endl;
    for(int i=0;i<solution.size();i++)
    {
        cout<<solution[i]<< " ";
     }cout<<endl;
     //cout<<"...........spider population initialization is taking place......................"<<endl;
    vector<Spider> population(10);
        Spider temp;
        temp=population[0];
        vector<int> t=max_degree();
        map<int,int> hash;
           for(int k=0;k<t.size();k++)
        hash[t[k]]=1;
    for(int i=0;i<10;i++)
    {   population[i].Position=t;
        for(int j=0;j<population[i].Position.size();j++)
        {   int done=0;
            while(!done){
                  //  cout<<"paglayo"<<endl;
                    double r=(double)(rand()/(RAND_MAX)+1.);
             int node=0+rand()%430;
            if(r>0.5)
                if(hash.find(node)==hash.end())
                {(population[i].Position)[j]=node;
                hash[(population[i].Position)[j]]=1;
                done=1;
                }
            }//cout<<"guoasd"<<endl;
        }
        population[i].tar.intensity=0;
       // cout<<"here";
        population[i].fitness=influence_set(population[i].Position);
       // cout<<"here 2"<<endl;
        population[i].cs=0;
        vector<int> mask(ss);
        vector<int> t(ss);
        population[i].dimension_mask=mask;
        population[i].prev_move=t;

     }/////////population init
   //  cout<<"...............population initialization done..............."<<endl;
 for(int p=0;p<100;p++){
        cout<<"iteration"<<p<<endl;
for(int i=0;i<10;i++)
{double max_vib=population[i].tar.intensity;
    for(int j=0;j<10;j++)
    {
        if(j!=i)
        {
             double res=source_vibration(population[j].fitness)*exp(-distance(population[i].Position,population[j].Position) / vibration_generation(population));
             //cout<<"res _vis"<<res<<endl;
             if(res>max_vib)
            {
                max_vib=res;

                population[i].pTarPostion=population[j].Position;
            }
        }
    }
    if(max_vib==population[i].tar.intensity)
    {population[i].cs++;
    }
    //cout<<"populatoon CS value"<<population[i].cs;}
    else
    {
        population[i].tar.intensity=max_vib;
        population[i].cs=0;
    }
    //cout<<"spider"<<i<<"population[i].tar.intensity"<<population[i].tar.intensity<<"target\t";
    for(int j=0;j<ss;j++)
    {
        cout<<population[i].Position[j]<<" ";   }cout<<endl;
        //cout<<" //////////////////////////////////////////////////////////////////////";
    double x=rand()/(RAND_MAX+1.0);
     if(x>pow(pc,population[i].cs))
    {
        int sum=0;
        for(int k=0;k<ss;k++)
        sum+=population[i].dimension_mask[k];
        if(sum==0)
        {
            int y= 1.0+rand()%ss;
            population[i].dimension_mask[y]=1;
        }
        else if(sum==ss){
              int y= 1.0+rand()%ss;
            population[i].dimension_mask[y]=0;
        }
        else{
            for(int k=0;k<ss;k++)
            {
            double y= rand()/(RAND_MAX+1.0);
            if(y>pm)
            population[i].dimension_mask[k]=0;
            else
             population[i].dimension_mask[k]=1;
            }
    }
    }
    vector<int> temp(ss);
     for(int k=0;k<ss;k++)
    {//cout<<k;
        if(!population[i].dimension_mask[k])
        {//cout<<k<<"inside if";
            temp[k]=population[i].pTarPostion[k];
        }
        else
        {//cout<<"inside else";
            int r= abs(rand()%10);//size of population
            //cout<<r;
            temp[k]=population[r].Position[k];
        }
    }
     double r=rand()/(RAND_MAX+1.0);
    map<int,int> hash2;
    for(int k=0;k<ss;k++){
        if(i==0)r=0;
        double R=(rand()%10)/10.0;
        int tantrum=population[i].Position[k];
         population[i].Position[k]=(abs((int)(population[i].Position[k]+(population[i].Position[k]-population[i].prev_move[k])*r+(temp[k]-population[i].Position[k])*R)%siz));
        int t=population[i].Position[k];
        if(hash2.find(population[i].Position[k])==hash2.end())
        {
            hash2[population[i].Position[k]]=1;
        }else
        {while(hash2.find(t)!=hash2.end())
        {
            t++;
        }population[i].Position[k]=t;}
        population[i].prev_move[k]=tantrum;
    }
    hash2.clear();
 }int index;
int max_fit=INT_MIN;
for(int i=0;i<10;i++)
{
    if(population[i].fitness>max_fit)
        {max_fit=population[i].fitness;
        index=i;
        }
}
cout<<"spider fitness";
cout<<population[index].fitness<<"\t ";
for(int i=0;i<ss;i++)
cout<<population[index].Position[i]<<" ";
cout<<endl;
cout<<endl;}
//
 t2=clock();
    float diff ((float)t2-(float)t1);
    cout<<diff<<endl;
    system ("pause");
//getch();
 }
