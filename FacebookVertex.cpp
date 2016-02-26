#include "FacebookVertex.h"

#include <vector>

FacebookVertex::FacebookVertex() : Vertex()
{
	name = "mynameisnobody";
}

FacebookVertex::FacebookVertex(string n) : Vertex()
{
	name = n;
}

string FacebookVertex::getName()
{
	return name;
}

void FacebookVertex::setName(string n)
{
	name = n;
}

string FacebookVertex::toString()
{
	string str;
	str = "V(";
	str.append(name);
	str.append(")");
	return str;
}

FacebookVertex* FacebookVertex::getVertex(string name, Graph* g)
{
	vector<Vertex*>* vertexList = g->vertices();
	for(int i = 0; i < vertexList->size(); i++)
	{
		if(name.compare(((FacebookVertex*)vertexList->at(i))->getName()) == 0)
		{
			return ((FacebookVertex*)vertexList->at(i));
		}
	}
	return NULL;
}

bool FacebookVertex::containsVertex(string name, Graph* g)
{
	if(getVertex(name, g) != NULL)
	{
		return true;
	}
	else
	{
		return false;
	}
}
