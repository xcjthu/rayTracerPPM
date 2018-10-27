#include "stdafx.h"
#include "Primitive.h"
#include <cmath>
#include <iostream>
#include <Eigen\Dense>


// using namespace Eigen;

Primitive::Primitive()
{
}


Primitive::~Primitive()
{
}


Color Texture::getColor(double u, double v) {
	int tmpi = int(bp.GetW() * abs(u));
	int tmpj = int(bp.GetH() * abs(v));
	// std::cout << tmpi << " " << tmpj << std::endl;
	return bp.GetColor(tmpj, tmpi);
}

double smallBox::Intersect(const Ray& ray) {
	double t;
	double dist = 1e20;
	t = (min.x - ray.o.x) / ray.dir.x;
	if (in(ray.o + ray.dir * t) && t < dist) dist = t;

	t = (min.y - ray.o.y) / ray.dir.y;
	if (in(ray.o + ray.dir * t) && t < dist) dist = t;

	t = (min.z - ray.o.z) / ray.dir.z;
	if (in(ray.o + ray.dir * t) && t < dist) dist = t;

	t = (max.x - ray.o.x) / ray.dir.x;
	if (in(ray.o + ray.dir * t) && t < dist) dist = t;

	t = (max.y - ray.o.y) / ray.dir.y;
	if (in(ray.o + ray.dir * t) && t < dist) dist = t;

	t = (max.z - ray.o.z) / ray.dir.z;
	if (in(ray.o + ray.dir * t) && t < dist) dist = t;

	return dist > 1e-4 ? dist : 1e20;
}

double Sphere::Intersect(const Ray& ray) {
	
	Vec op = o - ray.o;
	double t, b = op.dot(ray.dir);
	double det = b * b - op.dot(op) + radius * radius;
	if (det < 0) {
		return 1e20;
	}
	else {
		det = sqrt(det);
	}
	/*
	if (b - det > 1e-4) {
		t = b - det;
		crashPointLastTime = ray.o + ray.dir * t;
		return t;
	}
	else if (b + det > 1e-4) {
		t = b + det;
		crashPointLastTime = ray.o + ray.dir * t;
		return t;
	}
	return 1e20;
	*/
	return (t = b - det) > 1e-4 ? t : ((t = b + det)>1e-4 ? t : 1e20);
	
}

Color Sphere::getColor(Vec& pos) {
	if (texture == 0)
		return material.getColor();
	
	Vec det = pos - o;
	double tmpTheta = det.x / radius;
	double tmpAlpha = det.y / (radius * sin(acos(tmpTheta)));
	//std::cout << texture << std::endl;
	//std::cout << tmpAlpha << " " << tmpAlpha << std::endl;
	double tmpu = (tmpTheta + 1.) / 2;//(acos(tmpTheta) + 2 * 3.1416) / (4 * 3.1416);
	double tmpv = (tmpAlpha + 1.) / 2;// (acos(tmpAlpha) + 2 * 3.1416) / (4 * 3.1416);
	// std::cout << tmpu << " " << tmpv << std::endl;
	return texture->getColor(tmpu, tmpv);
}


double Plane::Intersect(const Ray& ray) {
	double d = nor.dot(ray.dir);
	if (d != 0) {
		double distTmp = -(nor.dot(ray.o) + D) / d;
		if (distTmp > 1e-4) {
			return distTmp;
		}
	}
	return 1e20;
}

Color Plane::getColor(Vec& pos) {
	if (texture == 0) return material.getColor();

	double tmpu;
	double tmpv;
	// tmpu = (pos - o)
	if (axis == 1) {
		tmpv = (pos.y * 0.01) - int(pos.y * 0.01);
		tmpu = (pos.z * 0.01) - int(pos.z * 0.01);
	}
	if (axis == 2) {
		tmpu = (pos.x * 0.01) - int(pos.x * 0.01);
		tmpv = (pos.z * 0.01) - int(pos.z * 0.01);
	}
	if (axis == 3) {
		tmpu = (pos.x * 0.01) - int(pos.x * 0.01);
		tmpv = (pos.y * 0.02) - int(pos.y * 0.02);
	}
	return texture->getColor(tmpu, tmpv);
}


Bezier::Bezier(int _m, int _n, std::vector<Vec> points, std::string textFile) : m(_m), n(_n) {

	for (int i = 0; i <= _m; ++i) {
		std::vector<Vec> tmp;
		control.push_back(tmp);
		for (int j = 0; j <= _n; ++j) {
			// std::cout << i * m + j << std::endl;
			control[i].push_back(points[i * (_m + 1) + j]);
			aabb.fit(points[i * (_m + 1) + j]);
		}
	}
	cn = new int[_n + 1];
	cn_1 = new int[_n];

	cm = new int[_m + 1];
	cm_1 = new int[_m];
	int tmpn = 1;
	
	cn[0] = cm[0] = cn_1[0] = cm_1[0] = 1;

	for (int i = 1; i <= _n; ++i) {
		cn[i] = cn[i - 1] * (_n - i + 1) / i;
		// std::cout << cn[i];
	}
	for (int i = 1; i < _n; ++i) {
		cn_1[i] = cn_1[i - 1] * (_n - i) / i;
		// std::cout << cn_1[i] << std::endl;
	}
	for (int i = 1; i <= _m; ++i) {
		cm[i] = cm[i - 1] * (_m - i + 1) / i;
	}
	for (int i = 1; i < _m; ++i) {
		cm_1[i] = cm_1[i - 1] * (_m - i) / i;
	}

	if (textFile == "None")
		texture = 0;
	else
		texture = new Texture(textFile);
}

Color Bezier::getColor(double u, double v) {
	
	if (texture == 0)
		return material.getColor();

	return texture->getColor(u, v);
}

Vec Bezier::getPoint(double u, double v) {
	
	Vec ans(0, 0, 0);
	for (int i = 0; i <= m; ++i) {
		for (int j = 0; j <= n; ++j) {
			ans = ans + control[i][j] * (cm[i] * pow(u, i) * pow(1 - u, m - i) * cn[j] * pow(v, j) * pow(1 - v, n - j));
		}
	}

	return ans;
}

void print(Vec a) { std::cout << a.x << " " << a.y << " " << a.z << std::endl; }
Vec Bezier::getDpDu(double u, double v) {
	Vec ansDpdu(0, 0, 0);
	for (int i = 1; i <= m; ++i) {
		double B1 = cm[i] * pow(u, i - 1) * pow(1 - u, m - i - 1) * (i - m * u);
		for (int j = 0; j <= n; ++j) {
			// double B1 = m * (cn_1[i - 1] * pow(u, i - 1) * pow(1 - u, n - i) - cm_1[i] * pow(u, i) * pow(1 - u, m - i - 1));
			double B2 = cn[j] * pow(v, j) * pow(1 - v, n - j);
			ansDpdu = ansDpdu + control[i][j] * (B1 * B2);
		}
	}
	for (int j = 0; j <= n; ++j) {
		// double B2 = cn[j] * pow(v, j) * pow(1 - v, n - j);
		ansDpdu = ansDpdu + control[0][j] * cn[j] * pow(v, j) * pow(1 - v, n - j) * (-m * pow(1 - u, m - 1));
	}
	return ansDpdu;
}

Vec Bezier::getDpDv(double u, double v) {
	Vec ans(0, 0, 0);
	for (int i = 0; i <= m; ++i) {
		double B1 = cm[i] * pow(u, i) * pow(1 - u, m - i);
		for (int j = 1; j <= n; ++j) {
			// double B2 = n * (cm_1[j - 1] * pow(v, j - 1) * pow(1 - v, m - j) - cn_1[j] * pow(v, j) * pow(1 - v, n - j - 1));
			double B2 = cn[j] * pow(v, j - 1) * pow(1 - v, n - j - 1) * (j - n * v);
			ans = ans + control[i][j] * (B1 * B2);
		}
		ans = ans + control[i][0] * cm[i] * pow(u, i) * pow(1 - u, m - i) * (-n * pow(1 - v, n - 1));
	}
	// for (int i = 0; i <= m; ++i) {
		// double B1 = cm[i] * pow(u, i) * pow(1 - u, m - i);
		// ans = ans + control[i][0] * cm[i] * pow(u, i) * pow(1 - u, m - i) * (-n * pow(1 - v, n - 1));
	// }
	return ans;
}

double Bezier::Intersect(const Ray& ray) {
	double dist = aabb.Intersect(ray);
	if (dist > 1e19) return 1e20;
	
	for (int index = 1; index < 20; ++index) {
		Vec initAns(dist, 0.05*index, 0.05*index);
		int num = 0;
		const int th = 20;

		while (num < th) {
			double t = initAns.x, u = initAns.y, v = initAns.z;
			// std::cout << "t: " << t << "u: " << u << "v: " << v << std::endl;
			
			Eigen::Vector3d p;
			p << t, u, v;
			Vec fx = ray.o + ray.dir * t - getPoint(u, v);
			num += 1;
			// std::cout << "|fx|£º" << fx.dot(fx) << std::endl;
			if (fx.dot(fx) < 1e-4) break;

			Eigen::Vector3d F;
			F << fx.x, fx.y, fx.z;
			Vec dpdu = getDpDu(u, v);
			Vec dpdv = getDpDv(u, v);
			/*
			for (int i = 1; i <= m; ++i) {
				for (int j = 0; j <= n; ++j) {
					double B1 = m * (cn_1[i - 1] * pow(u, i - 1) * pow(1 - u, n - i) - cm_1[i] * pow(u, i) * pow(1 - u, m - i - 1));
					double B2 = cn[j] * pow(v, j) * pow(1 - v, n - j);
					dpdu = dpdu + control[i][j] * (B1 * B2);
				}
			}*/
			/*
			for (int i = 0; i <= m; ++i) {
				for (int j = 1; j <= n; ++j) {
					double B1 = cm[i] * pow(u, i) * pow(1 - u, m - i);
					double B2 = n * (cm_1[j - 1] * pow(v, j - 1) * pow(1 - v, m - j) - cn_1[j] * pow(v, j) * pow(1 - v, n - j - 1));
					dpdv = dpdv + control[i][j] * (B1 * B2);
				}
			}*/

			Eigen::Matrix3d ma, mat;
			ma(0, 0) = ray.dir.x;
			ma(1, 0) = ray.dir.y;
			ma(2, 0) = ray.dir.z;

			ma(0, 1) = -dpdu.x;
			ma(1, 1) = -dpdu.y;
			ma(2, 1) = -dpdu.z;

			ma(0, 2) = -dpdv.x;
			ma(1, 2) = -dpdv.y;
			ma(2, 2) = -dpdv.z;

			mat = ma.inverse();
			// Eigen::Vector3d p1;
			p = p - mat * F;
			
			initAns.x = p(0);
			initAns.y = p(1);
			initAns.z = p(2);
		}

		double anst = initAns.x, ansu = initAns.y, ansv = initAns.z;
		if (num == th || ansu < 0 || ansu > 1 || ansv < 0 || ansv > 1) continue;

		// ans = t;
		Vec normal = getDpDu(ansu, ansv).crossProduct(getDpDv(ansu, ansv)).norm();
		if (normal.dot(ray.dir) < 0) apNor = normal;
		else apNor = normal * -1;
		lastCrashPointColor = getColor(ansu, ansv);
		return anst;
	}
	return 1e20;
	// return ans == dist ? 1e20 : ans;
}


Box::Box(std::string filepath, Material& ma, std::string textFile) {
	setMaterial(ma);
	FILE* f;
	fopen_s(&f, filepath.data(), "rt");
	int num;
	fscanf_s(f, "%d", &num);
	for (int i = 0; i < num; ++i) {
		int m, n;
		fscanf_s(f, "%d %d", &m, &n);
		//m += 1;
		//n += 1;
		std::vector<Vec> tmp;
		for (int j = 0; j < (m + 1) * (n + 1); ++j) {
			double tmpx, tmpy, tmpz;
			fscanf_s(f, "%lf %lf %lf", &tmpx, &tmpz, &tmpy);
			tmp.push_back(Vec(tmpx + 3, tmpy - 2, tmpz + 10) * 5 + 25);
			// Vec abcdef = Vec(tmpx, tmpy, tmpz) * 5 + 20;
			// std::cout << abcdef.x << " " << abcdef.y << " " << abcdef.z << std::endl;
		}
		Bezier* b = new Bezier(m, n, tmp, textFile);
		b->setMaterial(ma);
		//std::cout << i << std::endl;
		//print(b->getDpDu(0.1, 0.1));
		//print(b->getDpDv(0.1, 0.1));
		// print(b->getPoint(0.3, 0.3));
		// print(b->getPoint(0.5, 0.5));
		fit(b);
	}
	
}

/*
void Box::fit(const Vec& p) {
	if (p.x<min.x)min.x = p.x; // min
	if (p.y<min.y)min.y = p.y; // min
	if (p.z<min.z)min.z = p.z; // min
	max.x = p.x > max.x ? p.x : max.x;// MAX(p.x, max.x);
	max.y = p.y > max.y ? p.y : max.y;// MAX(p.y, max.y);
	max.z = p.z > max.z ? p.z : max.z;// MAX(p.z, max.z);
}
*/
void Box::fit(Bezier* thing) {
	pris.push_back(thing);
	for (int i = 0; i < thing->m; ++i) {
		for (int j = 0; j < thing->n; ++j) {
			// fit(thing->control[i][j]);
			aabb.fit(thing->control[i][j]);
		}
	}
}

double Box::Intersect(const Ray& ray) {
	double dist = aabb.Intersect(ray);
	if (dist > 1e19) return 1e20;

	double t = 1e20;
	for (int i = 0; i < pris.size(); ++i) {
		double tmp = pris[i]->Intersect(ray);
		if (tmp < t) {
			t = tmp;
			apNor = pris[i]->getNormal(Vec());
			lastCrashPointColor = pris[i]->lastCrashPointColor;
			//std::cout << t << std::endl;
		}
	}
	// return t;
	return t > 1e-4 ? t : 1e20;

}

