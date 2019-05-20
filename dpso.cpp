/*
Author -
IIT (BHU) Varanasi
Effective DPSO Implementation
*/
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <set>
#include <vector>
#define ve vector
using namespace std;
mt19937 Eng(chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count());
uniform_int_distribution<int> uid(0, 1000000000);
uniform_real_distribution<double> fid2(0, 1);
inline int rnd(int l, int r)
{
	// Generates random integer number in range [l,r]
	assert(l >= 0 && r >= l);
	return (uid(Eng) % (r - l + 1)) + l;
}
inline float frnd()
{
	// Generates a real number in the range [0.0,1.0]
	return fid2(Eng);
}
class dpso
{
	ve<int> Gbest;
	ve<double> prob;
	ve<ve<int> > v, X, Pbest;
	ve<ve<pair<int, float> >> adj;						 // adjacency list representation of graph containing (node,weight) pairs
	ve<long long> Aux_array1, Aux_array2, Aux_array3;		 // For intersection finding

	double w, c1, c2;
	long long marker;
	int n, sn, k, gmax;

public:
	dpso() :marker(1) {}

	// fitness value
	double LIE(ve<int> &);

	// Function to take input and allocate inital memory
	void input();

	// Function to find seed set
	ve<int> find_seed_set();

	// Function to initialize Pbest and Position X
	ve<ve<int> > degree_based_initialization();

	// Setting Aux_array[1,2] to find intersection with X[i]
	void vec_intersect(ve<int> &, ve<int> &);

	// modify passed vector with local best 
	void local_search(ve<int> &);

	void update_pos();

	void update_velocity();

	void update_gbest();

	void update_pbest();
};
int main()
{
	dpso d;
	d.input();
	return 0;
}

void dpso::input()
{
	/*
	Input format :
	n :		 Number of vertices in the graph
	sn:	     Number of Swarm Particles
	k :		 Size of seed set
	gmax:    Number of cycles
	c1,c2,w: constants required by algorithm
	m : number of edges

	next m lines describing edges (u,v,w)
	*/
	cin >> n >> sn >> k >> gmax >> c1 >> c2 >> w;
	Aux_array1.resize(n);
	Aux_array2.resize(n);
	Aux_array3.resize(n);
	prob.resize(n);
	v.resize(sn);

	// Each particle initialized to 0
	for (auto &x : v)
		x.resize(k);
	adj.resize(n);
	int m;
	cin >> m;
	for (int i = 0; i < m; i++)
	{
		int a, b;
		float c;
		cin >> a >> b >> c;
		--a;
		--b;
		adj[a].push_back({ b, c });
	}

	ve<int> seed_set = find_seed_set();
	cout << "Seed set is :\n";
	for (auto &x : seed_set)
		cout << x + 1 << " ";

	cout << "With fitness : " << LIE(seed_set) << endl;
	return;
}
ve<int> dpso::find_seed_set()
{
	X = degree_based_initialization();
	Pbest = degree_based_initialization();


	//step 2 choose Gbest
	update_gbest();

	// 3. begin cycling
	for (int g = 0; g < gmax; g++)
	{
		// update V velocity
		update_velocity();

		//update X position
		update_pos();

		//update Pbest
		update_pbest();

		//update_gbest
		auto copy_gbest = Gbest;
		double val1 = LIE(Gbest);

		update_gbest();

		// fine tune gbest using local search 
		local_search(Gbest);

		if (val1 > LIE(Gbest))
			Gbest = copy_gbest;
	}
	return Gbest;
}
void dpso::local_search(ve<int> &Xa)
{
	++marker;
	for (auto &x : Xa)
		Aux_array1[x] = marker;
	for (int i = 0; i < Xa.size(); i++)
	{
		bool flag = false;
		int node = Xa[i];
		ve<int> Neighbour;
		for (auto &x : adj[node])
		{
			if (Aux_array1[x.first] != marker)
				Neighbour.push_back(x.first);
		}
		if (Neighbour.empty())
			continue;
		random_shuffle(Neighbour.begin(), Neighbour.end());
		for (auto &y : Neighbour)
		{
			if (flag)
				break;
			int node2 = y;
			double val1 = LIE(Xa);
			Xa[i] = node2;
			double val2 = LIE(Xa);
			if (val1>val2)
			{
				Xa[i] = y;
			}
			else
			{
				flag = 1;
				Aux_array1[y] = marker - 1;
				Aux_array1[node2] = marker;
			}
		}
	}
	return;
}
double dpso::LIE(ve<int> &s)
{
	++marker;
	double ans = s.size();
	for (auto &x : s)
	{
		Aux_array1[x] = marker;
	}
	ve<int> hop1, hop2;
	for (auto &x : s)
	{
		for (auto &y : adj[x])
		{
			if (Aux_array1[y.first] == marker)
				continue;
			if (Aux_array2[y.first] != marker)
			{
				Aux_array2[y.first] = marker;
				prob[y.first] = (1 - y.second);
				hop1.push_back(y.first);
			}
			else
				prob[y.first] *= (1 - y.second);
		}
	}
	for (auto &x : hop1)
	{
		for (auto &y : adj[x])
		{
			if (Aux_array1[y.first] == marker || Aux_array2[y.first] == marker)
				continue;

			if (Aux_array3[y.first] != marker)
			{
				Aux_array3[y.first] = marker;
				prob[y.first] = (1 - y.second)*prob[x];
				hop2.push_back(y.first);
			}
			else
				prob[y.first] *= (1 - y.second)*prob[x];
		}
	}
	for (auto &x : hop1)
		ans += 1 - prob[x];
	for (auto &x : hop2)
		ans += 1 - prob[x];
	return ans;

}
void dpso::vec_intersect(ve<int> &v1, ve<int> &v2)
{
	++marker;
	for (auto &x : v1)
		Aux_array1[x] = marker;
	for (auto &x : v2)
		Aux_array2[x] = marker;
}
ve<ve<int> > dpso::degree_based_initialization()
{
	ve<int> idx(n);
	for (int i = 0; i < n; i++)
		idx[i] = i;
	sort(idx.begin(), idx.end(), [&](const int &a, const int &b) {return adj[a].size() > adj[b].size(); });
	ve<ve<int> > Z(sn);
	for (auto &x : Z)
		x.resize(k);
	for (int i = 0; i < sn; i++)
	{
		++marker;
		for (int j = 0; j < k; j++)
		{
			if (rnd(0, 1))
			{
				int repl = rnd(0, n - 1);
				while (Aux_array1[repl] == marker)
				{
					repl = rnd(0, n - 1);
				}

				Z[i][j] = repl;
				Aux_array1[repl] = marker;
			}
			else
			{
				int repl = idx[j];
				while (Aux_array1[repl] == marker)
				{
					repl = rnd(0, n - 1);
				}

				Z[i][j] = repl;
				Aux_array1[repl] = marker;
			}
		}
	}

	return Z;
}
void dpso::update_pos()
{
	for (int i = 0; i < sn; i++)
	{
		++marker;
		for (int j = 0; j < k; j++)
		{
			if (v[i][j] == 1)
			{
				int repl = rnd(0, n - 1);
				while (Aux_array1[repl] == marker)
				{
					repl = rnd(0, n - 1);
				}

				X[i][j] = repl;
				Aux_array1[repl] = marker;
			}
			else
			{
				int repl = X[i][j];
				while (Aux_array1[repl] == marker)
				{
					repl = rnd(0, n - 1);
				}

				X[i][j] = repl;
				Aux_array1[repl] = marker;
			}
		}
	}
}
void dpso::update_velocity()
{
	for (int i = 0; i < sn; i++)
	{
		vec_intersect(Pbest[i], Gbest);
		for (int j = 0; j < k; j++)
		{
			v[i][j] = (w * v[i][j] + c1 * frnd()*(Aux_array1[X[i][j]] != marker) +
				c2 * frnd()*(Aux_array2[X[i][j]] != marker)) >= 2;
		}
	}
}
void dpso::update_gbest()
{
	double val = 0;
	int id = 0;
	for (int i = 0; i < sn; i++)
	{
		double val2 = LIE(Pbest[i]);
		if (val2 > val)
			val = val2, id = i;
	}
	Gbest = Pbest[id];
}
void dpso::update_pbest()
{
	for (int i = 0; i < sn; i++)
	{
		for (int j = 0; j < k; j++)
		{
			if (rnd(0, 1) && adj[Pbest[i][j]].size())
			{
				int node = Pbest[i][j];
				Pbest[i][j] = adj[node][rnd(0, adj[node].size() - 1)].first;
			}
		}
	}
}