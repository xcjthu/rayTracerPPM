#pragma once

#include "Primitive.h"
#include "HashList.h"
#include "Bmp.h"
#include <cmath>

const double PI = 3.1415926;

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
	Bmp bp;
	Color* image;
	int num_photon = 1000000;
	
	bool Intersect(const Ray& ray, double& dist, int& id);
	void trace(const Ray &r, int dpt, int x, int y, double refrIndex, const Vec& adj);
	void trace0(const Ray &r, int dpt, bool m, const Vec &fl, const Vec &adj, int pixX, int pixY);
	void tracePass2(const Ray &r, int dpt, Vec& flux, double refrIndex);
	Vec randVec();
	void genp(Ray& r, Vec& flux);
	void render();
	void save();
	// tone mapping and gamma correction
	int toInt(double x) { return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5); }
	Engine(int _w = 1024, int _h = 768) : w(_w), h(_h), bp(_h, _w) { srand(19981102); scene = new Scene; }
	~Engine();
};

