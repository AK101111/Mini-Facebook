#ifndef GRAPH_H
#define GRAPH_H

#include<string>
#include<vector>
#include<map>

using namespace std;

class Vertex
{
	private:
		int vertex_number;
	public:
		Vertex();
		Vertex(int vertex_number);
		void setVertexNumber(int vertex_number);
		int getVertexNumber();
		string toString();
};

class Edge
{
	private:
		Vertex* source;
		Vertex* destination;
		bool directed;
		int weight;
//		int refs;
	public:
		Edge(Vertex* src, Vertex* dest);
		Edge(Vertex* src, Vertex* dest, bool isDir, int wt);
		Vertex* getSource();
		Vertex* getDestination();
		bool isDirected();
		int getWeight();
		string toString();
//		void setRefs(int r);
//		int getRefs();
};
	
class Graph
{
	private:
//		int num_vertices;
		int new_vertex_ID;
		vector<Vertex*>* vertexList;
		map<Vertex*, vector<Edge*>* >* adjacencyList;
		bool directed;
		bool weighted;
		Edge* removeAdjacencyListEdge(Vertex* v1, Vertex* v2);
	public:
		Graph();
		Graph(bool isDir, bool weighted);
		~Graph();
		void addVertex(Vertex* v);
		void addEdge(Vertex* src, Vertex* dest);
		void addEdge(Vertex* src, Vertex* dest, int wt);
		void addEdge(Vertex* src, Vertex* dest, int wt, bool isDir);
		void clear();
		void clearEdges();
		bool containsEdge(Edge* e);
		bool containsEdge(Vertex* v1, Vertex* v2);
		bool containsVertex(Vertex* v);
		int cost(vector<Vertex*> path);
		int degree(Vertex* v);
		Edge* edge(Vertex* v1, Vertex* v2);
		int edgeCount();
		vector<Edge*> edges();
		int edgeWeight(Vertex* v1, Vertex* v2);
		int inDegree(Vertex* v);
		bool isDirected();
		bool isEmpty();
		bool isReachable(Vertex* v1, Vertex* v2);		// DFS	
		bool isWeighted();
		vector<Vertex*>* minimumWeightPath(Vertex* v1, Vertex* v2);	// Dijkstra's algorithm
		vector<Vertex*> neighbors(Vertex* v);
		int outDegree(Vertex* v);
		Edge* removeEdge(Edge* e);
		Edge* removeEdge(Vertex* v1, Vertex* v2);
		void removeVertex(Vertex* v);
		vector<Vertex*>* shortestPath(Vertex* v1, Vertex* v2);	// BFS
		string toString();
		int vertexCount();
		vector<Vertex*>* vertices();
		map<Vertex*, vector<Edge*>* >* getAdjacencyList();
};

#endif
