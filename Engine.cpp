// #include "stdafx.h"
#include "Engine.h"

const double ALPHA = 0.7;

Engine::~Engine()
{
}


bool Engine::Intersect(const Ray& ray, double& dist, int& id) {
	int n = scene->numPrim;
	double d, inf = 1e20; 
	dist = inf;
	for (int i = 0; i<n; i++) {
		d = scene->getPrim(i)->Intersect(ray);
		if (d < dist) {
			dist = d;
			id = i;
		}
	}
	return dist < inf;
}


void Engine::trace(const Ray &r, int dpt, int x, int y, double refrIndex) {
	//�������� r: ����, dpt: ׷�����, m: �Ƿ���������֯����, f1:
	double dist;
	int id;

	++ dpt;
	if (!Intersect(r, dist, id) || dpt >= 20) return;

	Primitive* obj = scene->getPrim(id); //�ཻ������
	
	Vec ap = r.o + r.dir * dist; //�ཻ�Ľ���
	Vec n = obj->getNormal(ap).norm();

	Vec nl = n.dot(r.dir)<0 ? n : n*-1;

	// if (scene->getPrim(id)->getType() == "SPHERE") Sphere& obj
	
	if (obj->GetDiffuse() != 0) {
		//����һ��HitPoint
		hashL.addNode(HPoint(ap, n, r.dir, x, y));
		// HPoint(const Vec _pos, const Vec _nor, const Vec _rD, int _x, int _y):pos(_pos),nor(_nor),rayDir(_rD), x(_x), y(_y){}

	}
	if (obj->GetReflection() != 0) {
		//����������·��
		trace(Ray(ap, r.dir - n * 2.0 * n.dot(r.dir)), x, y, dpt, refrIndex);
	}
	if (obj->GetRefraction() != 0) {
		//����������·��
		double refr = obj->GetRefraction();
		if (refr > 0) {
			double objRefr = obj->GetRefrIndex();
			
			Ray lr(ap, r.dir - n*2.0*n.dot(r.dir));
			bool into = (n.dot(nl)>0.0);
			double ddn = r.dir.dot(nl), cos2t;
			double nnt = into ? refrIndex / objRefr : objRefr / refrIndex;
			// total internal reflection
			if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) //ȫ����ֻ�з������
				return trace(lr, dpt, x, y, refrIndex);

			Vec td = (r.dir * nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
			// double a = objRefr - refrIndex, b = objRefr + refrIndex;
			// double R0 = a*a / (b*b);
			// double c = 1 - (into ? -ddn : td.dot(n));
			// double Re = R0 + (1 - R0)*c*c*c*c*c;
			Ray rr(ap, td);
			
			//������ߺ�������߶���
			trace(lr, dpt, x, y, refrIndex);
			trace(rr, dpt, x, y, objRefr);
		}
	}
}


void Engine::tracePass2(const Ray& r, int dpt, int x, int y, Vec& flux, double refrIndex){
	double dist;
	int id;
	
	++ dpt;
	if (!Intersect(r, dist, id) || dpt >= 20) return;
	
	Primitive* obj = scene->getPrim(id); //�ཻ������
	
	Vec ap = r.o + r.dir * dist; //�ཻ�Ľ���
	Vec n = obj->getNormal(ap).norm();
	Color f = obj->getColor();
	
	Vec nl = n.dot(r.dir)<0 ? n : n*-1;
	
	if (obj->GetDiffuse() != 0){
		Vec hh = (ap - hashL.hpBox.min) * hashL.hash_s;
		int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
		{
			std::vector<HPoint*> hp = hashL.hashGrid[hashL.hash(ix, iy, iz)];
			for(int tmp = 0; tmp < hp.size(); ++ tmp){
				Vec v = hp[tmp]->pos - ap;
				if (hp[tmp]->nor.dot(n) > 1e-3 && v.dot(v) <= hp[tmp]->radius2){
					double g = (hp[tmp]->N * ALPHA + ALPHA)/(hp[tmp]->N * ALPHA + 1.0);
					hp[tmp]->radius2 = hp[tmp]->radius2 * g;
					++ hp[tmp]->N;
					// hp[tmp]->flux = (hp[tmp]->flux + hp[tmp]->f.mul(fl)*(1./PI))*g??
				}
			}
			// tracePass2(<#const Ray &r#>, <#int dpt#>, <#int x#>, <#int y#>, <#Vec &flux#>, <#double refrIndex#>)??
		}
	}
	if (obj->GetReflection() != 0){
		tracePass2(Ray(ap, r.dir - n * 2.0 * n.dot(r.dir)), x, y, dpt, f.mul(flux),refrIndex);
	}
	if (obj->GetRefraction() != 0){
		// ??
	}
	
	
}

void Scene::initScene() {
	prim = new Primitive*[20];
	/*
	prim[0] = new Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), DIFF),//Left
	prim[1] = new Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), DIFF),//Right
	prim[2] = new Sphere(1e5, Vec(50, 40.8, 1e5), Vec(.75, .75, .75), DIFF),//Back
	prim[3] = new Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), DIFF),//Front
	prim[4] = new Sphere(1e5, Vec(50, 1e5, 81.6), Vec(.75, .75, .75), DIFF),//Bottomm
	prim[5] = new Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), DIFF),//Top
	prim[6] = new Sphere(16.5, Vec(27, 16.5, 47), Vec(1, 1, 1)*.999, SPEC),//Mirror
	prim[7] = new Sphere(16.5, Vec(73, 16.5, 88), Vec(1, 1, 1)*.999, REFR),//Glass
	prim[8] = new Sphere(8.5, Vec(50, 8.5, 60), Vec(1, 1, 1)*.999, DIFF)
	*/
}
