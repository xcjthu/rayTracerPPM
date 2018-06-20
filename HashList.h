#pragma once

#include "Vec.h"

#define MAX(x, y) ((x > y) ? x : y)
struct HPoint {
	
	Vec pos, nor, rayDir;
	int  brdfIndex;
	int x, y; //ÏñËØÎ»ÖÃ
	Color wgt;
	double radius;
	int N; //photon count
	Color flux;

};


class HashNode {
public:
	HashNode(HPoint* i, HashNode* _next) { hp = i; next = _next; };
	~HashNode();

	HPoint* hp;
	HashNode* next;
	// List* listAdd(HPoint*i);
};


class AABB {//°üÎ§ºÐ
public:
	Vec min, max; // axis aligned bounding box
	void fit(const Vec &p)
	{
		if (p.x<min.x)min.x = p.x; // min
		if (p.y<min.y)min.y = p.y; // min
		if (p.z<min.z)min.z = p.z; // min
		max.x = MAX(p.x, max.x);
		max.y = MAX(p.y, max.y);
		max.z = MAX(p.z, max.z);
	}
	void reset() {
		min = Vec(1e20, 1e20, 1e20);
		max = Vec(-1e20, -1e20, -1e20);
	}
};

class HashList {
public:
	HashNode* head;
	HashNode* tail;

	HashList() { head = tail = 0; }
	
	void addNode(HPoint* hp){}
};

/*
struct List { HPoint *id; List *next; };
List* ListAdd(HPoint *i, List* h) {
	List* p = new List;
	p->id = i;
	p->next = h;
	return p;
}
*/
