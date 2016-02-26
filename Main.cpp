#include<iostream>
#include<string>
#include<stdlib.h>
#include<fstream>
#include<cstring>

#include "Graph.h"
#include "FacebookVertex.h"

using namespace std;

void populateGraph(Graph* g, string filename);
void printFacebookGraph(Graph* g);
void erdos(Graph* g);
void reachable(Graph* g);
void shortestPath(Graph* g);

const int ACQUAINTANCE = 3;
const int FRIEND = 2;
const int CLOSE_FRIEND = 1;


int main(int argc, char **argv)
{

	cout << "COL 106 Assignment 6 : Facebook" << endl;
	cout << "===============================" << endl;

	if(argc != 2)
	{
		cout << endl << "Usage : ./a.out <absolute-path-to-graphfile>" << endl;
		exit(1);
	}

        cout << endl << "input filename : " << argv[1] << endl;

	Graph* g = new Graph();
	populateGraph(g, argv[1]);

	char resp = 'y';
	while(resp == 'y')
	{
		cout << "\nenter 1 to compute erdos number of a user : " << endl;
		cout << "enter 2 to check if a user can reach another through a series of friends : " << endl;
		cout << "enter 3 to compute the shortest path between two users : " << endl;
		char c;
		cin >> c;
		switch(c)
		{
			case '1' : {
					erdos(g);
					break;
				}
			case '2' : {
					reachable(g);
					break;
				}
			case '3' : {
					shortestPath(g);
					break;
				}
		}
		cout << "do you wish to continue? ";
		cin >> resp;
	}

	delete(g);

	return 0;
}

void erdos(Graph* g)
{
	string username;
	cout << "enter username : ";
	cin >> username;

	FacebookVertex* v = FacebookVertex::getVertex(username, g);
	if(v == NULL)
	{
		cout << "no such user : " << username << endl;
		return;
	}

	FacebookVertex* erdos_v = FacebookVertex::getVertex("Erdos", g);
	if(erdos_v == NULL)
	{
		cout << "no such user : Erdos" << endl;
		return;
	}

	vector<Vertex*>* shortestPath = g->shortestPath(erdos_v, v);

	if(shortestPath == NULL)
	{
		cout << "no path exists between Erdos and " << username << endl;
	}
	else
	{
		cout << "erdos number = " << shortestPath->size()-1 << endl;
		delete(shortestPath);
	}
}

void reachable(Graph* g)
{
	string username1, username2;
	cout << "enter username 1 : ";
	cin >> username1;
	cout << "enter username 2 : ";
	cin >> username2;

	FacebookVertex* v1 = FacebookVertex::getVertex(username1, g);
	if(v1 == NULL)
	{
		cout << "no such user : " << username1 << endl;
		return;
	}

	FacebookVertex* v2 = FacebookVertex::getVertex(username2, g);
	if(v2 == NULL)
	{
		cout << "no such user : " << username2 << endl;
		return;
	}

	cout << "is " << username1 << " reachable from " << username2 << " : ";
	if(g->isReachable(v1, v2) == true)
	{
		cout << "true" << endl;
	}
	else
	{
		cout << "false" << endl;
	}
}

void shortestPath(Graph* g)
{
	string username1, username2;
	cout << "enter username 1 : ";
	cin >> username1;
	cout << "enter username 2 : ";
	cin >> username2;

	FacebookVertex* v1 = FacebookVertex::getVertex(username1, g);
	if(v1 == NULL)
	{
		cout << "no such user : " << username1 << endl;
		return;
	}

	FacebookVertex* v2 = FacebookVertex::getVertex(username2, g);
	if(v2 == NULL)
	{
		cout << "no such user : " << username2 << endl;
		return;
	}

	vector<Vertex*>* path = g->minimumWeightPath(v1, v2);
	
	if(path == NULL)
	{
		cout << "no path exists" << endl;
	}
	else
	{
		cout << "minimum weight path : " << endl;
		for(int i = 0; i < path->size(); i++)
		{
			cout << ((FacebookVertex*)path->at(i))->toString() << " ";
		}
		cout << endl << "path weight = " << g->cost(*path) << endl;
		delete(path);
	}
}

void populateGraph(Graph* g, string filename)
{
	fstream f;
	char *buf;
	string sbuf;

	f.open(filename.data(), ios::in);
	if(f == NULL)
	{
		string err_msg;
		err_msg += "opening ";
		err_msg += filename;
		err_msg += " : FAILED!";
		cout << err_msg << endl;
		exit(1);
	}

	char * pch;
	getline(f, sbuf);
	while(!f.eof())
	{
		if(sbuf.size() > 0)
		{
			buf = new char[sbuf.size()+1];
			strcpy(buf, sbuf.c_str());

			//friend 1
			FacebookVertex* f1 = NULL;
			string f1name;
			pch = strtok (buf," \t\n");
			f1name = pch;
			if(FacebookVertex::containsVertex(f1name, g) == false)
			{
				f1 = new FacebookVertex(f1name);
				g->addVertex(f1);
			}
			f1 = FacebookVertex::getVertex(f1name, g);

			pch = strtok (NULL," \t\n");
			if(pch != NULL)
			{
				//friend 2
				FacebookVertex* f2 = NULL;
				string f2name;
				f2name = pch;
				if(FacebookVertex::containsVertex(f2name, g) == false)
				{
					f2 = new FacebookVertex(f2name);
					g->addVertex(f2);
				}
				f2 = FacebookVertex::getVertex(f2name, g);
				
				//edge
				pch = strtok (NULL, " \t\n");
				string edgeWeight = pch;
				if(edgeWeight.compare("Acquaintance") == 0)
				{
					g->addEdge(f1, f2, ACQUAINTANCE);
				}
				else if(edgeWeight.compare("Friend") == 0)
				{
					g->addEdge(f1, f2, FRIEND);
				}
				else if(edgeWeight.compare("CloseFriend") == 0)
				{
					g->addEdge(f1, f2, CLOSE_FRIEND);
				}
				else
				{
					cout << "unknown edge weight" << endl;
					cout << "exiting.." << endl;exit(1);
				}
			}

			delete(buf);
		}
		getline(f, sbuf);
	}

	printFacebookGraph(g);
}

void printFacebookGraph(Graph* g)
{
	vector<Vertex*>* vertexList = g->vertices();
	map<Vertex*, vector<Edge*>* >* adjacencyList = g->getAdjacencyList();

	string str;
	str.append("\n----------------------------------------------------------------------------\n\nGraph:\n\n");
	for(int i = 0; i < vertexList->size(); i++)
	{
		Vertex* v = vertexList->at(i);
		str.append(((FacebookVertex*)v)->toString());
		str.append(" connected to\n");
		for(int j = 0; j < adjacencyList->at(v)->size(); j++)
		{
			if(v == adjacencyList->at(v)->at(j)->getSource())
			{
				str.append(((FacebookVertex*)adjacencyList->at(v)->at(j)->getDestination())->toString());
			}
			else if(v == adjacencyList->at(v)->at(j)->getDestination())
			{
				str.append(((FacebookVertex*)adjacencyList->at(v)->at(j)->getSource())->toString());
			}
			else
			{
				cout << "Function : string Graph::toString()\n";
				cout << "unrelated edge : " << adjacencyList->at(v)->at(j)->toString() << " in adjacency list of " << ((FacebookVertex*)v)->toString() << endl;
				cout << "exiting.." << endl;exit(1);
			}
			str.append(" [W=");
			char buf[32];
			sprintf(buf, "%d", adjacencyList->at(v)->at(j)->getWeight());
			str.append(buf);
			str.append("]    ");
		}
		str.append("\n\n");
	}
	str.append("\n----------------------------------------------------------------------------\n\n");
	cout << str;
}

void test(Graph* g)
{
/*	cout << "hello world\n";

	Graph* g = new Graph();
	Vertex* v1 = new Vertex();
	g->addVertex(v1);
	Vertex* v2 = new Vertex();
	g->addVertex(v2);
	Vertex* v3 = new Vertex();
	g->addVertex(v3);
	Vertex* v4 = new Vertex();
	g->addVertex(v4);
	g->addEdge(v1, v2, 1, true);
	g->addEdge(v2, v3, 2);
	g->addEdge(v3, v4, 3);
	cout << g->toString();

//	g->removeVertex(v2);
//	cout << g->toString();

	g->removeEdge(v3, v1);
	cout << g->toString();

	cout << "indegree of V(1) = " << g->inDegree(v2) << endl;
	cout << "outdegree of V(1) = " << g->outDegree(v2) << endl;
	cout << "degree of V(1) = " << g->degree(v2) << endl;
	cout << "indegree of V(0) = " << g->inDegree(v1) << endl;
	cout << "outdegree of V(0) = " << g->outDegree(v1) << endl;
	cout << "degree of V(0) = " << g->degree(v1) << endl;

	Vertex* v5 = new Vertex();
	cout << "containsEdge(V(0), V(1)) = " << g->containsEdge(v1, v2) << endl;
	cout << "containsEdge(V(1), V(0)) = " << g->containsEdge(v2, v1) << endl;
	cout << "containsEdge(V(0), V(2)) = " << g->containsEdge(v1, v3) << endl;
	cout << "containsEdge(V(2), V(0)) = " << g->containsEdge(v3, v1) << endl;
	cout << "containsEdge(V(1), V(2)) = " << g->containsEdge(v2, v3) << endl;
	cout << "containsEdge(V(2), V(1)) = " << g->containsEdge(v3, v2) << endl;
	cout << "containsEdge(V(0), V(4)) = " << g->containsEdge(v1, v5) << endl;
//	cout << "containsEdge(V(4), V(0)) = " << g->containsEdge(v5, v1) << endl;
//	cout << "containsEdge(V(4), V(4)) = " << g->containsEdge(v5, v5) << endl;

	vector<Vertex*> t;
	t.push_back(v1);
	t.push_back(v2);
	t.push_back(v3);
	t.push_back(v4);
	cout << "cost = " << g->cost(t) << endl;

	g->removeEdge(v2, v3);
	cout << g->toString() << endl;

	g->clearEdges();
	cout << g->toString();
	g->removeVertex(v3);
	cout << "V(2) present? : " << g->containsVertex(v3) <<endl;
	cout << g->toString();
	g->addVertex(v3);
	cout << "V(2) present? : " << g->containsVertex(v3) <<endl;
	cout << g->toString();
	g->clear();
	cout << g->toString();*/
}
