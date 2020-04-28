//verification.cpp
//part of the PESPlib library
//written by Marc Goerigk
//Institute for Numerical and Applied Mathematics
//University of Goettingen, Germany


#include<set>
#include<vector>
#include<list>
#include<map>
#include<iostream>
#include<fstream>
#include<cassert>
#include<cstdlib>
#include<queue>

using namespace std;

struct Neighbor
{
	int n;
	int l;
	int id;
	int order;
	
	Neighbor(int _n, int _l, int _id=0, int _order=1) : n(_n), l(_l), id(_id), order(_order) {}
	
	bool operator< (const Neighbor& other) const
	{
	  return (n < other.n);
	}
};

struct Edge
{
	int from;
	int to;
	int weight;
	int lower;
	int upper;
	int id;

	Edge(int f, int t, int w, int l, int u, int i) : from(f), to(t), weight(w), lower(l), upper(u), id(i){}	

	bool operator< (const Edge& other) const
	{
		return (id < other.id);
	}
};

typedef set<Edge> Graph;

class Kruskal
{
public:
	void solve(Graph g);
	vector<list<Neighbor> > get_tree();
	Graph get_edgeset();
private:
	Graph t;
};

class EAN
{
public:
	void check();
	void read_activities(char* file);
	void read_timetable(char* file);
	void write(char* file);
private:
	Graph g;
	int num_vertices;
	int num_edges;
	map<int,int> lower_map;
	map<int,int> timetable;
};

void Kruskal::solve(Graph g)
{
	list<set<int> > nodes;
	list<set<int> >::iterator l,r;
	
	for (Graph::iterator g_it = g.begin(); g_it!=g.end(); ++g_it)
	{
		l = nodes.end();
		r = nodes.end();
		
		for (list<set<int> >::iterator it=nodes.begin(); it!=nodes.end(); ++it)
		{
			if (it->find(g_it->from) != it->end())
				l = it;
			if (it->find(g_it->to) != it->end())
				r = it;
		}
		
		if (!(l==r && l!=nodes.end()))
		{
			t.insert(*g_it);
			if (l!=nodes.end() && r!=nodes.end())
			{
				l->insert(r->begin(),r->end());
				nodes.erase(r);
			}
			else if(l!=nodes.end() && r==nodes.end())
			{
				l->insert(g_it->from);
				l->insert(g_it->to);
			}
			else if(r!=nodes.end() && l==nodes.end())
			{
				r->insert(g_it->from);
				r->insert(g_it->to);
			}
			else
			{
				set<int> temp;
				temp.insert(g_it->from);
				temp.insert(g_it->to);
				nodes.push_back(temp);
			}
		}
	}
}

vector<list<Neighbor> > Kruskal::get_tree()
{
      int num_vertices=0;
      for (Graph::iterator it=t.begin(); it!=t.end(); ++it)
      {
	  num_vertices = std::max(num_vertices, std::max(it->from, it->to));
      }
      
      vector<list<Neighbor> > sol(num_vertices);
      for (Graph::iterator it=t.begin(); it!=t.end(); ++it)
      {
	  sol[it->from-1].push_back(Neighbor(it->to-1, it->lower, it->id, 1));
	  sol[it->to-1].push_back(Neighbor(it->from-1, it->lower, it->id, -1));
      }
      
      return sol;
}

Graph Kruskal::get_edgeset()
{
    return t;
}

void EAN::check()
{ 	
	bool feas = true;
	
	for (set<Edge>::iterator it=g.begin(); it!=g.end(); ++it)
		if (it->lower > timetable[it->id] || it->upper < timetable[it->id])
		{
			cout<<"Infeasibility detected: Edge "<<it->id<<"\n";
			feas = false;
		}
	
	Kruskal kruskal;
 	kruskal.solve(g);
	
	vector<list<Neighbor> > tree = kruskal.get_tree();
	Graph intree = kruskal.get_edgeset();
	Graph outtree;
	
	for (Graph::iterator it=g.begin(); it!=g.end(); ++it)
	{  
	    if (intree.count(*it)==0)
	      outtree.insert(*it);
	}
	
	cout<<"Calculating cycles...   ";
	int loopcounter = 0;
	
	list<list<int> > cycles;
	for (Graph::iterator outedge=outtree.begin(); outedge != outtree.end(); ++outedge)
	{
	    cycles.push_back(list<int>());
	    (cycles.back()).push_back(-outedge->id);
	    
	    queue<int> q;
	    set<int> marks;
	    map<int,list<int> > path;
	    
	    q.push(outedge->from);
	    marks.insert(outedge->from);
	    list<int> emptylist;
	    path[outedge->from] = emptylist;
	    
	    bool stilltogo = true;
	    while (q.size() > 0 && stilltogo)
	    {
			int v = q.front();
			q.pop();
			if (v == outedge->to)
			    stilltogo = false;
			else
			for (list<Neighbor>::iterator it=tree[v-1].begin(); it!=tree[v-1].end(); ++it)
			{
			    int w = it->n+1;
			    if (marks.count(w)==0)
			    {
				marks.insert(w);
				q.push(w);
				path[w] = path[v];
				if (it->order == 1)
				  path[w].push_back(it->id);
				else
				  path[w].push_back(-it->id);
			    }
			}
		
	    }
	    
	    for (list<int>::iterator it=path[outedge->to].begin(); it!=path[outedge->to].end(); ++it)
			(cycles.back()).push_back(*it);
	    
	    ++loopcounter;
	    int perc = 100 - 100*(outtree.size() - loopcounter) / outtree.size();
	    if ( perc < 10 )
			cout<<"\b\b"<< perc<<"%"<<flush;
		else if (perc<100)
			cout<<"\b\b\b"<< perc<<"%"<<flush;
		else
			cout<<"\b\b\b\b"<< perc<<"%"<<flush;
	    
	}
	
	cout<<"\n";
	
	for (list<list<int> >::iterator it=cycles.begin(); it!=cycles.end(); ++it)
	{	
		int cycletension = 0;
	    for (list<int>::iterator it2=it->begin(); it2!=it->end(); ++it2)
	    {
			if (*it2>0)
				cycletension += timetable[*it2];
			else
				cycletension -= timetable[-*it2];
	    }
	    cycletension%=60;
	    if (cycletension != 0)
	    {
			cout<<"Infeasibility detected in cycle ";
			for (list<int>::iterator it2=it->begin(); it2!=it->end(); ++it2)
				cout<<*it2<<" ";
		    cout<<"\n";
		    feas = false;
		}
	}
	
	if (feas)
		cout<<"Timetable is feasible.\n";
		
	int objective = 0;
	for (Graph::iterator it=g.begin(); it!=g.end(); ++it)
	{  
	    objective += (timetable[it->id] - lower_map[it->id]) * it->weight;
	}
	
	cout<<"Objective is "<<objective<<"\n";
	
}

void EAN::read_activities(char* file)
{	
	num_vertices = 0;
	num_edges = 0;
	
    ifstream activities(file);
    assert(activities != 0);
    string line;
    int num_heads = 0;
    while (!activities.eof())
    {
        getline(activities,line);
        if (line == "" || (line.c_str())[0]<48 || (line.c_str())[0]>57)
            continue;

	int id, tail, head, min, max, weight;
	
        //activity-id
	id = atoi(line.c_str());
        size_t pos = line.find(";");
        line=line.substr(pos+1);
	  
        //tail
        tail=atoi(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);

        //head
        head=atoi(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);

        //lower
        min=atoi(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);
	
        //upper
        max=atoi(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);	

        //passengers
        weight=atoi(line.c_str());
        
	  g.insert(Edge(tail,head,weight,min,max,id));
	  num_vertices = std::max(std::max(tail,head),num_vertices);
	  num_edges = std::max(num_edges, id);
	  lower_map[id] = min;
    }
    
    activities.close();	
	
}

void EAN::read_timetable(char* file)
{
    ifstream timfile(file);
    assert(timfile != 0);
    string line;
    int num_heads = 0;
    while (!timfile.eof())
    {
        getline(timfile,line);
        if (line == "" || (line.c_str())[0]<48 || (line.c_str())[0]>57)
            continue;

		int id, duration;
	
        //activity-id
		id = atoi(line.c_str());
        size_t pos = line.find(";");
        line=line.substr(pos+1);
	  
        //duration
        duration=atoi(line.c_str());
        pos = line.find(";");
        line=line.substr(pos+1);

     timetable[id]=duration;
     
    }
    
    timfile.close();	
}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		cout<<"Usage: ./verification INSTANCEFILE SOLUTIONFILE\n";
		exit(0);
	}
	
	EAN ean;
	ean.read_activities(argv[1]);
	ean.read_timetable(argv[2]);
	ean.check();
	
	return 0;
}
