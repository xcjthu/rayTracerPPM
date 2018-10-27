#pragma once

#include "Vec.h"
#include <string>
#include <vector>
#include "Bmp.h"
#include <iostream>


class Material {
public:
	Color color;
	double diff;
	double refl;
	double refr;
	double refrIndex;
	// double specular;

	//double drefl;

public:
	void SetColor(Color& c) { color = c; }
	Color getColor() { return color; }
	void SetDiffuse(double a_Diff) { diff = a_Diff; }
	void SetReflection(double a_Refl) { refl = a_Refl; }
	void SetRefraction(double a_Refr) { refr = a_Refr; }
	void SetRefrIndex(double r_index) { refrIndex = r_index; }
	// void SetSpecular(double a) { specular = a; }
	// double GetSpecular() { return specular; }
	double GetSpecular() { return 1 - diff; }
	double GetDiffuse() { return diff; }
	double GetReflection() { return refl; }
	double GetRefraction() { return refr; }
	double GetRefrIndex() { return refrIndex; }
	Material(const Color c = Color(0, 0, 0), const double _diff = 0, double _refl = 0, double _refr = 0) :color(c), diff(_diff), refl(_refl), refr(_refr) { refrIndex = 0; };
	~Material() {};
};


class Texture {
public:
	Bmp bp;
	Color getColor(double u, double v); // u,v参数
	Texture(std::string inputImg) { bp.Input(inputImg); }
};

struct smallBox {
	Vec min, max;
	smallBox() { reset(); }
	void fit(const Vec& p) {
		if (p.x<min.x)min.x = p.x; // min
		if (p.y<min.y)min.y = p.y; // min
		if (p.z<min.z)min.z = p.z; // min
		max.x = p.x > max.x ? p.x : max.x;// MAX(p.x, max.x);
		max.y = p.y > max.y ? p.y : max.y;// MAX(p.y, max.y);
		max.z = p.z > max.z ? p.z : max.z;// MAX(p.z, max.z);
	}
	void reset() { min = Vec(1., 1., 1.) * 1e20; max = Vec(1., 1., 1.) * -1e20; }
	bool in(const Vec& p) { return (p.x - min.x >= -1e-4) && (p.x - max.x <= 1e-4) && (p.y - min.y >= -1e-4) && (p.y - max.y <= 1e-4) && (p.z - min.z >= -1e-4) && (p.z - max.z <= 1e-4); }
	double Intersect(const Ray& ray);
};

class Primitive
{
public:
	Material material;
	char* name;
	bool mLight;

	Primitive();
	~Primitive();

	Material* getMaterial() { return &material; }
	void setMaterial(Material& mat) { material = mat; }

	virtual double Intersect(const Ray& a_Ray) = 0;
	virtual Vec getNormal(Vec& pos) = 0;
	virtual void Light(bool a_light) { mLight = a_light; }
	virtual Color getColor(Vec& pos) { return material.getColor(); }
	virtual std::string getType() { return "Primitive"; };
	double GetSpecular() { return material.GetSpecular(); }
	double GetDiffuse() { return material.GetDiffuse(); }
	double GetReflection() { return material.GetReflection(); }
	double GetRefraction() { return material.GetRefraction(); }
	double GetRefrIndex() { return material.GetRefrIndex(); }
	bool isLight() { return mLight; }
	char* getName() { return name; }
	void setName(char* aName) {
		name = new char[strlen(aName) + 1];
		for (int i = 0; i < strlen(aName); ++i) {
			name[i] = aName[i];
		}
		name[strlen(aName)] = '\0';
	}
};


class Sphere :public Primitive {
public:
	double radius;
	Vec o;

	Texture* texture;


	Sphere(const Vec& _o, double r, Material& m, std::string textfile = "None") :o(_o), radius(r) { 
		setMaterial(m);  
		if (textfile == "None")
			texture = 0;
		else
			texture = new Texture(textfile);
	}
	// void setTexture(Texture* tmpt) { texture = tmpt; }

	~Sphere() {}
	// Vec crashPointLastTime;
	Color getColor(Vec& pos);

	double Intersect(const Ray& a_Ray);
	Vec getNormal(Vec& pos) { return (pos - o).norm(); }
	std::string getType() { return "SPHERE"; }
};


class Plane :public Primitive {
public:
	Vec nor;
	double D;

	Texture* texture;
	
	Plane(const Vec& n, double d, Material& m, std::string textFile = "None") :nor(n), D(d) { 
		if (n.x != 0) axis = 1;
		else if (n.y != 0) axis = 2;
		else axis = 3;
		setMaterial(m); 
		if (textFile == "None")
			texture = 0;
		else
			texture = new Texture(textFile);
	}
	int axis = 0;
	Color getColor(Vec& pos);
	double Intersect(const Ray& a_Ray);
	Vec getNormal(Vec& pos) { return nor; }
	std::string getType() { return "PLANE"; }
};


class Bezier : public Primitive {
public:
	int m, n;//曲面有w*h个控制点
	Texture* texture;//纹理
	smallBox aabb;//包围盒
	std::vector<std::vector<Vec>> control;//控制点
	//提前计算好的一些参数
	int* cn;
	int* cn_1;
	int* cm;
	int* cm_1;

	Vec apNor;//上一次求交的交点
	Vec getPoint(double u, double v);//根据参数求曲面上的点
	Vec getNormal(Vec& pos) { return apNor; }//求法向量
	Vec getDpDu(double u, double v);//求导
	Vec getDpDv(double u, double v);
	Color getColor(double u, double v);//获取点对应的纹理的颜色
	Color lastCrashPointColor;//上一次求交的交点的颜色
	double Intersect(const Ray& a_ray);//求交，并且返回交点
	Bezier(int _m, int _n, std::vector<Vec> points, std::string textFile = "None");
};


class Box:public Primitive {
public:
	std::vector<Bezier*> pris;
	Vec apNor;

	smallBox aabb;
	Color lastCrashPointColor;
	// Vec min, max; // axis aligned bounding box
	void fit(Bezier* thing);
	double Intersect(const Ray& ray);
	
	Vec getNormal(Vec& pos) { return apNor; }
	Color getColor(Vec& pos) { return lastCrashPointColor; }
	// void fit(const Vec &p);
	
	// void reset() { min = Vec(1e20, 1e20, 1e20); max = Vec(-1e20, -1e20, -1e20); }

	Box(std::string filePath, Material& ma, std::string textFile = "None");
	Box(Material& ma) { setMaterial(ma); }
};


