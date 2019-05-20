#include <functional>
#include <iostream>
#include <algorithm>
#include <vector>
#include <ctime>
#include <set>
#include <map>
#define ve vector
#pragma warning(disable:4996)
#define pi pair<int,ve<double>>
#define pdi pair<double,int>
using namespace std;
int readInt();
class best_Effort
{
	/**
		n-nodes count
		m-edges count
		t-topic count
		Y-current query distribution
		adj-adjacency list of graph
		rev- reverse adjacency list of a graph
		edges : set of edges of graph
		vis : marker to check whether a node has been visited or not1
		state : state of a node
		seed : set of nodes currently selected as seeds
		L,H : Priority queues contaning the upper_bound of influence margin of nodes
	*/
	int n, m, t, k;
	ve<double> Y;
	ve<ve<pi> > adj;
	map<int, map<int, ve<double>> > edges;

	multiset<pair<double, int>, greater<pair<double, int> > > L, H;
	ve<int> vis;
	ve<int> state;       // 0-initial,2-exact
	set<int> seed;
	ve<ve<pi> > rev;
	int version;
	ve<int> run()
	{
		L.clear();
		H.clear();
		seed.clear();
		fill(state.begin(),state.end(),0);
		for (int i = 0; i<n; i++)
			L.insert({ est_inf(i),i });
		clock_t c = clock();
		//clock_t e=c+CLOCKS_PER_SEC*1;
		int N=n;
		ve<int> id(n);
		for(auto &x:id)
			x=&x-&id[0];
		sort(id.begin(),id.end(),[&](const int &a,const int &b){return adj[a].size()<adj[b].size();});
		for (int i = 0; i<k; i++)
		{
			seed.insert(id[i]);
		}
		
		printf("\nExpected number of activated nodes : %.3lf\n",expected_activates());
		return ve<int>(seed.begin(), seed.end());
	}
	double pp(int u, int v) // u-influencing v;
	{
		double ans = 0;
		if (edges[u][v].size() == 0)
			return 0;
		for (int i = 0; i<t; i++)
		{
			ans += Y[i] * edges[u][v][i];
		}
		return ans;
	}
	double activate_uv(int u, int v)
	{
		++version;
		set<pair<double, int>,  greater<pair<double, int> >> U;
		U.insert({ 0,u });
		while (!U.empty())
		{
			int beg = U.begin()->second;
			double val = U.begin()->first;
			// v is reached return the cost incurred to reach v
			if (beg == v)
				return U.begin()->first;
			vis[beg] = version;
			for (auto &x : adj[beg])
			{
				if (vis[x.first] != version)
				{
					vis[x.first] = version;
					U.insert({ val + pp(beg,x.first),x.first });
				}
			}
		}
		// v is unreachable from u hence return 0
		return 0;
	}
	double ap(int u)
	{
		++version;
		return apt(u);
	}
	double apt(int u)
	{
		if (seed.find(u) != seed.end())
			return 1;
		double api = 0, z = 1;
		vis[u] = version;
		for (auto &x : rev[u])
		{
			if (vis[x.first] != version)
			{
				vis[x.first] = version;
				z *= (1 - apt(x.first)*pp(x.first, u));
			}
		}
		return 1 - z;
	}
	double expected_activates()
	{
		double val = 0;
		for (int i = 0; i<n; i++)
			val += ap(i);
		return val;
	}
	double margin(int u)
	{
		double d1 = expected_activates();
		seed.insert(u);
		d1 = expected_activates() - d1;
		seed.erase(u);
		return d1;
	}
public:

	void input();
	void query();
	double est_inf(int u)
	{
		++version;
		return bound_calc(u);
	}
	double bound_calc(int u)
	{
		vis[u] = version;
		double bv = 0;
		for (auto &x : adj[u])
		{
			if (vis[x.first] != version)
			{
				bv += *max_element(x.second.begin(), x.second.end());
				bv += bound_calc(x.first);
			}
		}
		return bv;
	}
	void insert_candidates()
	{
		while (!L.empty() && (H.empty() || *H.begin() <= *L.begin()))
		{
			H.insert(*L.begin());
			L.erase(L.begin());
		}
	}
};
int main()
{
	/*
		Input :
			n-graph nodes count
			m-edges count
			t-topic count
			m triplets (u,v,w) : Edges u to v with influence w { w is a vector }
			Q-Number of TIM queries
			then Q TIM of form [K,Y] in two lines each
		Sample input :
			3 2 2
			1 2 0.2 0.1
			2 3 0.1 0.3
			2
			1
			0.2 0.1
			2
			0.5 0.8
		Sample output :
			Expected number of activated nodes : 1.052
			Seeds :
			1
			Expected number of activated nodes : 2.290
			Seeds :
			1 2
	*/
	best_Effort B;
	B.input();
	int K=readInt();
	for(int i=0;i<K;i++)
		B.query();
	return 0;
}
int readInt()
{
	int n;
	scanf("%d", &n);
	return n;
}
void best_Effort::input()
{
	
	n = readInt();
	m = readInt();
	t = readInt();
	state.resize(n);
	adj.resize(n);
	vis.resize(n);
	rev.resize(n);
	fill(vis.begin(), vis.begin() + n, 0);
	version = 0;
	for (int i = 0; i<m; i++)
	{
		int x = readInt() - 1, y = readInt() - 1;
		int c1 = adj[x].size();
		int c2 = rev[y].size();
		adj[x].push_back({ y,ve<double>() });
		rev[y].push_back({ x,ve<double> ()});
		for (int j = 0; j<t; j++)
		{
			double d;
			scanf("%lf", &d);
			adj[x][c1].second.push_back(d);
			rev[y][c2].second.push_back(d);
		}
		edges[x][y] = adj[x][c1].second;
	}

}
void best_Effort::query()
{
	k = readInt();
	Y.resize(t);
	for (int i = 0; i<t; i++)
		scanf("%lf", &Y[i]);
	// create ordering in decreasing order
	/*for (int i = 0; i<n; i++)
	{
		sort(adj[i].begin(), adj[i].end(), [&](const pi &a, const pi &b)
		{
			return pp(i, a.first)>pp(i, b.first);
		});
		sort(rev[i].begin(), rev[i].end(), [&](const pi &a, const pi &b)
		{
			return pp(a.first, i)>pp(b.first, i);
		});
	}*/
	ve<int> vi(run());
	cout << "Seeds :\n";
	for (int i = 0; i<vi.size(); i++)
	{
		cout << vi[i]+1 << " ";
	}
	cout << endl;
}
