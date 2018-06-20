#pragma once

#include "Primitive.h"


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

	bool Intersect(const Ray& ray, double& dist, int& id);
	void trace(const Ray &r, int dpt, double refrIndex);
	Engine(int _w = 1024, int _h = 768) : w(_w), h(_h) {}
	~Engine();
};

