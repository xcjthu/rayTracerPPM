#pragma once

#include "Primitive.h"
#include "HashList.h"


class Scene {
public:
	Primitive** prim;
	int numPrim;

	void initScene();
	Primitive* getPrim(int index) { return prim[index]; }
	Scene() :prim(0), numPrim(0) { initScene(); }
	~Scene() { delete prim; }
};

class Engine
{
public:
	Scene* scene;
	int w, h;
	HashList hashL;
	
	bool Intersect(const Ray& ray, double& dist, int& id);
	void trace(const Ray &r, int dpt, int x, int y, double refrIndex);
	void tracePass2(const Ray &r, int dpt, int x, int y, Vec& flux, double refrIndex);
	Vec randVec(){return Vec(rand(), rand(), rand()).norm();}
	Engine(int _w = 1024, int _h = 768) : w(_w), h(_h) {srand(19981102);}
	~Engine();
};

