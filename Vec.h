#pragma once

#include <cmath>


class Vec{
public:
	double x, y, z; // vector: position, also color (r,g,b)
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }
	Vec(const Vec& a) { x = a.x; y = a.y; z = a.z; }
	inline Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	inline Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	inline Vec operator+(double b) const { return Vec(x + b, y + b, z + b); }
	inline Vec operator-(double b) const { return Vec(x - b, y - b, z - b); }
	inline Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
	inline Vec mul(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	inline Vec norm() { return (*this) * (1.0 / sqrt(x*x + y*y + z*z)); }
	inline double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
	Vec operator%(Vec&b) { return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
};

typedef Vec Color;
/*
class Color : public Vec {

};
*/

class Ray {
public:
	Ray(const Vec _o, const Vec _dir) : o(_o), dir(_dir) {}
	const Vec& getOri() { return o; }
	const Vec& getDir() { return dir; }
	Ray reflect(const Vec& n) { return Ray(x, r.d - n*2.0*n.dot(r.d); }
	Vec o, dir;
};

