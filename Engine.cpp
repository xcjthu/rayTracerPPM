#include "stdafx.h"
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
	//参数意义 r: 光线, dpt: 追踪深度, m: 是否是用来组织光子, f1:
	double dist;
	int id;

	++ dpt;
	if (!Intersect(r, dist, id) || dpt >= 20) return;

	Primitive* obj = scene->getPrim(id); //相交的物体
	
	Vec ap = r.o + r.dir * dist; //相交的交点
	Vec n = obj->getNormal(ap).norm();

	Vec nl = n.dot(r.dir)<0 ? n : n*-1;

	// if (scene->getPrim(id)->getType() == "SPHERE") Sphere& obj
	
	if (obj->GetDiffuse() != 0) {
		//建立一个HitPoint
		hashL.addNode(HPoint(ap, n, r.dir, x, y));
		// HPoint(const Vec _pos, const Vec _nor, const Vec _rD, int _x, int _y):pos(_pos),nor(_nor),rayDir(_rD), x(_x), y(_y){}

	}
	if (obj->GetReflection() != 0) {
		//按照正常光路走
		trace(Ray(ap, r.dir - n * 2.0 * n.dot(r.dir)), x, y, dpt, refrIndex);
	}
	if (obj->GetRefraction() != 0) {
		//按照正常光路走
		double refr = obj->GetRefraction();
		if (refr > 0) {
			double objRefr = obj->GetRefrIndex();
			
			Ray lr(ap, r.dir - n*2.0*n.dot(r.dir));
			bool into = (n.dot(nl)>0.0);
			double ddn = r.dir.dot(nl), cos2t;
			double nnt = into ? refrIndex / objRefr : objRefr / refrIndex;
			// total internal reflection
			if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) //全反射只有反射光线
				return trace(lr, dpt, x, y, refrIndex);

			Vec td = (r.dir * nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
			// double a = objRefr - refrIndex, b = objRefr + refrIndex;
			// double R0 = a*a / (b*b);
			// double c = 1 - (into ? -ddn : td.dot(n));
			// double Re = R0 + (1 - R0)*c*c*c*c*c;
			Ray rr(ap, td);
			
			//反射光线和折射光线都有
			trace(lr, dpt, x, y, refrIndex);
			trace(rr, dpt, x, y, objRefr);
		}
	}
}

Vec Engine::randVec() {
	double cost1 = 2 * PI * (double(rand())) / RAND_MAX;
	double cost2 = 2 * PI * (double(rand())) / RAND_MAX;
	return Vec(cos(cost1), sin(cost1) * cos(cost2), sin(cost1) * sin(cost2));
}

void Engine::genp(Ray& ray, Vec& flux) {
	flux = Vec(2500, 2500, 2500) * (PI * 4.0);
	ray.dir = randVec();
	ray.o = Vec(50, 50, 50); //光源

}

void Engine::tracePass2(const Ray& r, int dpt, Vec& flux, double refrIndex){
	double dist;
	int id;
	
	++ dpt;
	if (!Intersect(r, dist, id) || dpt >= 20) return;
	
	Primitive* obj = scene->getPrim(id); //相交的物体
	
	Vec ap = r.o + r.dir * dist; //相交的交点
	Vec n = obj->getNormal(ap).norm();
	Color f = obj->getColor();

	double p = f.x>f.y&&f.x>f.z ? f.x : f.y>f.z ? f.y : f.z;//吸收的概率
	
	Vec nl = n.dot(r.dir)<0 ? n : n*-1;

	double diffPar = obj->GetDiffuse();
	double reflPar = obj->GetReflection();
	double refrPar = obj->GetRefraction();

	double randNum = (diffPar + reflPar + refrPar + p) * double(rand()) / RAND_MAX;
	int choice;
	if (randNum < diffPar) choice = 0;
	else if (randNum < diffPar + reflPar) choice = 1;
	else choice = 2;

	
	// if (obj->GetDiffuse() != 0){
	if(choice == 0){
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
					hp[tmp]->flux = (hp[tmp]->flux + hp[tmp]->wgt.mul(flux)*(1. / PI))*g;
				}
			}
			tracePass2(Ray(ap, randVec()), dpt, f.mul(flux)*(1. / p), refrIndex);
		}
	}
	// if (obj->GetReflection() != 0){
	if (choice == 1){
		tracePass2(Ray(ap, r.dir - n * 2.0 * n.dot(r.dir)), dpt, f.mul(flux),refrIndex);
	}
	// if (obj->GetRefraction() != 0){
	if (choice == 2){
		// ??
		Ray lr(ap, r.dir - n*2.0*n.dot(r.dir));
		bool into = (n.dot(nl)>0.0);
		double nc = refrIndex, nt = obj->GetRefrIndex(), nnt = into ? nc / nt : nt / nc, ddn = r.dir.dot(nl), cos2t;

		// total internal reflection
		if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) return tracePass2(lr, dpt, flux, refrIndex);

		Vec td = (r.dir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
		double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : td.dot(n));
		double Re = R0 + (1 - R0)*c*c*c*c*c, P = Re; Ray rr(ap, td); 
		
		// photon ray (pick one via Russian roulette)
		tracePass2(rr, dpt, flux, nt);
	}	
}

void Engine::render() {
	int w = 1024, h = 768;
	// trace eye rays and store measurement points
	Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	Vec cx = Vec(w*.5135 / h), cy = (cx % cam.dir).norm()*.5135, *c = new Vec[w*h], vw;
	for (int y = 0; y<h; y++) {
		fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y / (h - 1));
		for (int x = 0; x<w; x++) {
			// pixel_index = x + y * w;
			Vec d = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5) + cam.dir;
			trace(Ray(cam.o + d * 140, d.norm()), 0, x, y, 1.0);
		}
	}

	hashL.build_hash_grid(w, h);

	int num_photon = 1000;

	vw = Vec(1, 1, 1);
	// #pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i<num_photon; i++) {
		double p = 100.*(i + 1) / num_photon;
		fprintf(stderr, "\rPhotonPass %5.2f%%", p);
		int m = 1000 * i;
		Ray r;
		Vec f;
		for (int j = 0; j<1000; j++) {
			genp(r, f);
			tracePass2(r, 0, f, 1.0);
		}
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
