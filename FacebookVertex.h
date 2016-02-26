#ifndef FACEBOOK_VERTEX_H
#define FACEBOOK_VERTEX_H

#include<string>
#include<iostream>
#include"Graph.h"

using namespace std;

class FacebookVertex : public Vertex
{
	private:
		string name;
	public:
		FacebookVertex();
		FacebookVertex(string name);
		string getName();
		void setName(string name);
		string toString();
		static FacebookVertex* getVertex(string name, Graph* g);
		static bool containsVertex(string name, Graph* g);
};

#endif
