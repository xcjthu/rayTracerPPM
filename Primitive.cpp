#include "stdafx.h"
#include "Primitive.h"


Primitive::Primitive()
{
}


Primitive::~Primitive()
{
}


double Sphere::Intersect(const Ray& ray) {
	/*
	Vec op = o - ray.o;
	double t, b = op.dot(ray.dir);
	double det = b * b - op.dot(op) + radius * radius;
	if (det < 0) {
		return 1e20;
	}
	else {
		det = sqrt(det);
	}
	return (t = b - det) > 1e-4 ? t : ((t = b + det)>1e-4 ? t : 1e20);
	*/
	
	Vec v = ray.o;
	v = v - o;
	Vec dir = ray.dir;

	double b = -v.dot(ray.dir);
	double det = b * b - v.dot(v) + radius * radius;
	if (det > 0) {
		det = sqrtf(det);
		double i1 = b - det;
		double i2 = b + det;
		if (i2 > 0) {
			if (i1 < 0) {
				return i2;
			}
			else {
				return i1;
			}
		}
	}
	return 1e20;
	
}


double Plane::Intersect(const Ray& ray) {
	double d = nor.dot(ray.dir);
	if (d != 0) {
		double distTmp = -(nor.dot(ray.o) + D) / d;
		if (distTmp > 0) {
			return distTmp;
		}
	}
	return 1e20;
}