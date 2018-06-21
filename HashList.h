#pragma once

#include "Vec.h"
#include <vector>

#define MAX(x, y) ((x > y) ? x : y)
struct HPoint {
	
	Vec pos, nor, rayDir;
	int x, y; //ÏñËØÎ»ÖÃ
	double radius2;
	int N; //photon count
	Color wgt;
	int  brdfIndex;
	Color flux;
	
	HPoint(const Vec _pos, const Vec _nor, const Vec _rD, int _x, int _y):pos(_pos),nor(_nor),rayDir(_rD), x(_x), y(_y){}
};

/*
class HashNode {
public:
	HashNode(HPoint* i, HashNode* _next) { hp = i; next = _next; };
	~HashNode();

	HPoint* hp;
	HashNode* next;
	// List* listAdd(HPoint*i);
};
*/

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
	std::vector<HPoint*> hList;
	AABB hpBox;
	double hash_s;
	int num_hash;
	std::vector<HPoint*>* hashGrid;
	
	HashList(){}
	
	/*
	void addNode(HPoint* hp){
	}*/
	void addNode(HPoint hp){ hList.push_back(&hp); }
	void build_hash_grid(int w, int h);
	int hash(const int ix, const int iy, const int iz){
		return (unsigned int)((ix*73856093)^(iy*19349663)^(iz*83492791)) % num_hash;
	}
};



