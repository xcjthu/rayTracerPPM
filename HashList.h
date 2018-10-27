#pragma once

#include "Vec.h"
#include <vector>

#define MAX(x, y) ((x > y) ? x : y)
struct HPoint {
	
	Vec pos, nor, rayDir;
	int x, y; //像素位置
	double radius2;
	int N; //photon count
	Color wgt;
	int  brdfIndex;
	Color flux;
	//Vec f;
	
	HPoint(const Vec _pos, const Vec _nor, const Vec _rD, const Vec _f, int _x, int _y):pos(_pos),nor(_nor),rayDir(_rD), x(_x), y(_y), wgt(_f){}
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

class AABB {//非常简单的包围盒，仅仅用于计算普通的大小。
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
	std::vector<HPoint*> hList;//所有的HitPoint
	AABB hpBox;//包围盒
	double hash_s;
	int num_hash;//Hash表的大小
	std::vector<HPoint*>* hashGrid;//每个网格中对应的HitPoint链表
	void clear();//清空
	HashList() { hashGrid = 0; }
	void addNode(HPoint* hp){ hList.push_back(hp); }//添加HitPoint
	void build_hash_grid(int w, int h);//在找到所有的HitPoint后建立Hash网格
	int hash(const int ix, const int iy, const int iz){//Hash值
		return (unsigned int)((ix*73856093)^(iy*19349663)^(iz*83492791)) % num_hash;
	}
};



