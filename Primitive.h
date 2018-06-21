#pragma once

#include "Vec.h"
#include <string>
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
	Material(const Color c = Color(0, 0, 0), const double _diff = 0, double _refl = 0, double _refr = 0):color(c), diff(_diff), refl(_refl), refr(_refr){};
	~Material() {};
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
	virtual Color getColor() { return material.getColor(); }
	virtual std::string getType() = 0;
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

	Sphere(const Vec& _o, double r, Material& m) :o(_o), radius(r) { setMaterial(m); }
	~Sphere() {}

	double Intersect(const Ray& a_Ray);
	Vec getNormal(Vec& pos) { return (pos - o).norm(); }
	std::string getType() { return "SPHERE"; }
};


class Plane :public Primitive {
public:
	Vec nor;
	double D;
	
	Plane(const Vec& n, double d, Material& m) :nor(n), D(d) { setMaterial(m); }
	double Intersect(const Ray& a_Ray);
	Vec getNormal(Vec& pos) { return nor; }
	std::string getType() { return "PLANE"; }
};
