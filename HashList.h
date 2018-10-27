#pragma once

#include "Vec.h"
#include <vector>

#define MAX(x, y) ((x > y) ? x : y)
struct HPoint {
	
	Vec pos, nor, rayDir;
	int x, y; //����λ��
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

class AABB {//�ǳ��򵥵İ�Χ�У��������ڼ�����ͨ�Ĵ�С��
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
	std::vector<HPoint*> hList;//���е�HitPoint
	AABB hpBox;//��Χ��
	double hash_s;
	int num_hash;//Hash��Ĵ�С
	std::vector<HPoint*>* hashGrid;//ÿ�������ж�Ӧ��HitPoint����
	void clear();//���
	HashList() { hashGrid = 0; }
	void addNode(HPoint* hp){ hList.push_back(hp); }//���HitPoint
	void build_hash_grid(int w, int h);//���ҵ����е�HitPoint����Hash����
	int hash(const int ix, const int iy, const int iz){//Hashֵ
		return (unsigned int)((ix*73856093)^(iy*19349663)^(iz*83492791)) % num_hash;
	}
};



