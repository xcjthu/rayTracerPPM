#include "stdafx.h"
#include "Engine.h"
#include <iostream>

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

double hal(int d) {
	return ((double)rand()) / RAND_MAX;
}

void Engine::trace(const Ray &r, int dpt, int x, int y, double refrIndex, const Vec& adj) {
	//参数意义 r: 光线, dpt: 追踪深度
	double dist;
	int id;

	++ dpt;
	if (!Intersect(r, dist, id) || dpt >= 20) return;
	
	// std::cout << "intersect\n" << std::endl;

	Primitive* obj = scene->getPrim(id); //相交的物体
	Color f = obj->getColor();

	Vec ap = r.o + r.dir * dist; //相交的交点
	Vec n = obj->getNormal(ap).norm();

	Vec nl = n.dot(r.dir)<0 ? n : n*-1;

	// if (scene->getPrim(id)->getType() == "SPHERE") Sphere& obj
	
	if (obj->GetDiffuse() > 0) {
		//建立一个HitPoint
		HPoint* tmp = new HPoint(ap, n, r.dir, f.mul(adj), x, y);
		hashL.addNode(tmp);
		// result = true;
		
	}
	if (obj->GetReflection() > 0) {
		//按照正常光路走
		trace(Ray(ap, r.dir - n * (2.0 * n.dot(r.dir))), x, y, dpt, refrIndex, f.mul(adj));
		// if (result) std::cout << "nothing" << std::endl;
	}
	if (obj->GetRefraction() > 0) {
		//按照正常光路走
		double refr = obj->GetRefraction();
		if (refr > 0) {
			double objRefr = obj->GetRefrIndex();
			
			Ray lr(ap, r.dir - n*2.0*n.dot(r.dir));
			bool into = (n.dot(nl)>0.0);
			// double ddn = r.dir.dot(nl), cos2t;
			// double nnt = into ? refrIndex / objRefr : objRefr / refrIndex;
			double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.dir.dot(nl), cos2t;
			// total internal reflection
			if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) //全反射只有反射光线
				return trace(lr, dpt, x, y,  refrIndex, adj);

			Vec td = (r.dir * nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
			// double a = objRefr - refrIndex, b = objRefr + refrIndex;
			double a = nt - nc, b = nt + nc;
			double R0 = a*a / (b*b);
			double c = 1 - (into ? -ddn : td.dot(n));
			double Re = R0 + (1 - R0)*c*c*c*c*c;
			Ray rr(ap, td);
			Vec fa = f.mul(adj);
			
			//反射光线和折射光线都有
			trace(lr, dpt, x, y, refrIndex, fa*Re);
			trace(rr, dpt, x, y, objRefr, fa*(1.0-Re));

			// if (result) std::cout << "abcdefg" << std::endl;
		}
	}
}



void Engine::trace0(const Ray &r, int dpt, bool m, const Vec &fl, const Vec &adj, int pixX, int pixY)
{
	double t;
	int id;

	dpt++;
	if (!Intersect(r, t, id) || (dpt >= 20))return;

	int d3 = dpt * 3;
	Primitive* obj = scene->getPrim(id);
	Vec x = r.o + r.dir*t;//交点
	// n = (x - obj->p).norm(), //交点处圆的法向量
	Vec n = obj->getNormal(x);
	Vec f = obj->getColor();//？？
	Vec nl = n.dot(r.dir)<0 ? n : n*-1; //调整方向后的法向量
	double p = f.x>f.y&&f.x>f.z ? f.x : f.y>f.z ? f.y : f.z;//f中最大值，RGB最大值？
	// if (p == 0) p = 0.001;
	double random = ((double)rand()) / RAND_MAX;


	if (obj->GetDiffuse() != 0) {//漫散射
		// Lambertian

		// use QMC to sample the next direction
		double r1 = 2.*PI*hal(d3 - 1), r2 = hal(d3 + 0); //随机数
		double r2s = sqrt(r2);

		Vec w = nl;
		Vec u = ((fabs(w.x)>.1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w%u, d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();


		if (m) {
			// eye ray
			// store the measurment point
			HPoint* hp = new HPoint(x, n, r.dir, f.mul(adj), pixX, pixY);
			hashL.addNode(hp);
		}
		else
		{
			// photon ray
			// find neighboring measurement points and accumulate flux via progressive density estimation
			Vec hh = (x - hashL.hpBox.min) * hashL.hash_s;
			int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
			{
				// List* hp = hash_grid[hash(ix, iy, iz)];
				std::vector<HPoint*> hp = hashL.hashGrid[hashL.hash(ix, iy, iz)];
				for (int index = 0; index < hp.size(); ++index) {
					Vec v = hp[index]->pos - x;
					if ((hp[index]->nor.dot(n) > 1e-3) && (v.dot(v) <= hp[index]->radius2)) {
						double g = (hp[index]->N * ALPHA + ALPHA) / (hp[index]->N*ALPHA + 1.0);
						hp[index]->radius2 = hp[index]->radius2 * g;
						hp[index]->N++;
						hp[index]->flux = (hp[index]->flux + hp[index]->wgt.mul(fl)*(1.0/ PI))*g;
						// std::cout << "flux: (" << hp[index]->flux.x << "," << hp[index]->flux.y << "," << hp[index]->flux.z << "," << fl.x << ")" << std::endl;
						
					}
				}
			}
			if (hal(d3 + 1)<p) trace0(Ray(x, d), dpt, m, f.mul(fl)*(1./ p), adj, pixX, pixY);
		}

	}
	if (obj->GetReflection() > 0) {
		// mirror
		trace0(Ray(x, r.dir - n*2.0*n.dot(r.dir)), dpt, m, f.mul(fl), f.mul(adj), pixX, pixY);

	}
	if (obj->GetRefraction() > 0) {
		// glass
		Ray lr(x, r.dir - n*2.0*n.dot(r.dir));
		bool into = (n.dot(nl)>0.0);
		double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.dir.dot(nl), cos2t;

		// total internal reflection
		if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) return trace0(lr, dpt, m, fl, adj, pixX, pixY);

		Vec td = (r.dir*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
		double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : td.dot(n));
		double Re = R0 + (1 - R0)*c*c*c*c*c, P = Re; Ray rr(x, td); Vec fa = f.mul(adj);
		if (m) {
			// eye ray (trace both rays)
			trace0(lr, dpt, m, fl, fa*Re, pixX, pixY);
			trace0(rr, dpt, m, fl, fa*(1.0 - Re), pixX, pixY);
		}
		else {
			// photon ray (pick one via Russian roulette)
			
			(hal(d3 - 1)<P) ? trace0(lr, dpt, m, fl, fa, pixX, pixY) : trace0(rr, dpt, m, fl, fa, pixX, pixY);
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
	ray.o = Vec(50, 60, 85); //光源

}

void Engine::tracePass2(const Ray& r, int dpt, Vec& flux, double refrIndex){
	double dist;
	int id;
	
	++ dpt;
	// std::cout << "dpt:" << dpt << std::endl;
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
		tracePass2(Ray(ap, r.dir - n * 2.0 * n.dot(r.dir)), dpt, f.mul(flux), refrIndex);
	}
	// if (obj->GetRefraction() != 0){
	if (choice == 2){
		// 
		Ray lr(ap, r.dir - n*2.0*n.dot(r.dir));
		bool into = (n.dot(nl)>0.0);
		// double nc = refrIndex, nt = obj->GetRefrIndex(), nnt = into ? nc / nt : nt / nc, ddn = r.dir.dot(nl), cos2t;
		double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.dir.dot(nl);
		double cos2t;
		// if (nt > 0)std::cout << nt << std::endl;
		// total internal reflection
		if ((cos2t = 1 - nnt * nnt * (1 - ddn*ddn))<0) return tracePass2(lr, dpt, flux, refrIndex);

		Vec td = (r.dir * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
		// double a = nt - nc;
		// double b = nt + nc;
		// double R0 = a*a / (b*b);
		// double c = 1 - (into ? -ddn : td.dot(n));
		// double Re = R0 + (1 - R0)*c*c*c*c*c; 
		Ray rr(ap, td);
		
		// photon ray (pick one via Russian roulette)
		tracePass2(rr, dpt, flux, nt);
	}	
}


double toDouble(double x) {
	return pow(1 - exp(-x), 1 / 2.2);
}

void Engine::save() {

	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			image[j * w + i].x = image[j * w + i].y = image[j * w + i].z = 0;
		}
	}

	for (int i = 0; i < hashL.hList.size(); ++i) {
		int tmpx = hashL.hList[i]->x;
		int tmpy = hashL.hList[i]->y;
		image[tmpy * w + tmpx] = image[tmpy * w + tmpx] + hashL.hList[i]->flux * (1.0 / (PI * hashL.hList[i]->radius2 * num_photon * 1000.0));

	}

	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			bp.SetColor(i, j, Color(toDouble(image[i * w + j].x), toDouble(image[i * w + j].y), toDouble(image[i * w + j].z)));
		}
	}

	bp.Output("D:/learn/graphics/image.bmp");
}

void Engine::render() {
	// int w = 1024, h = 768;

	image = new Color[w * h + 1];

	// trace eye rays and store measurement points
	Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	Vec cx = Vec(w*.5135 / h), cy = (cx % cam.dir).norm()*.5135, *c = new Vec[w*h], vw;
	for (int y = 0; y<h; y++) {
		fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y / (h - 1));
		for (int x = 0; x<w; x++) {
			// pixel_index = x + y * w;
			bool tmpb = false;
			Vec d = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5) + cam.dir;
			// trace(Ray(cam.o + d * 140, d.norm()), 0, x, y, 1.0, Vec(1, 1, 1));
			trace0(Ray(cam.o + d * 140, d.norm()), 0, true, Vec(), Vec(1, 1, 1), x, y);
			// std::cout << hashL.hList.size() << ":" << hashL.hList[hashL.hList.size() - 1]->flux.x << std::endl;
		}
	}

	hashL.build_hash_grid(w, h);
	std::cout << "\nsize:" << hashL.hList.size() << std::endl;
	
	/*
	for (int i = 0; i < hashL.hList.size(); ++i) {
		std::cout << "flux: (" << hashL.hList[i]->flux.x << "," << hashL.hList[i]->flux.y << "," << hashL.hList[i]->flux.z << ")" << std::endl;
		std::cout << "x:" << hashL.hList[i]->x << " y:" << hashL.hList[i]->y << std::endl;
	}*/

	// int num_photon = 1000;
	

	vw = Vec(1, 1, 1);
	// #pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i<num_photon; i++) {
		double p = 100.*(i + 1) / num_photon;
		// std::cout << i << std::endl;
		fprintf(stderr, "\rPhotonPass %5.2f%%", p);
		//std::cout << "abcdefg" << std::endl;
		int m = 1000 * i;
		Ray r;
		Vec f;
		for (int j = 0; j<1000; j++) {
			
			genp(r, f);
			trace0(r, 0, 0>1, f, vw, 0, 0);
			//tracePass2(r, 0, f, 1.0);
		}
		if (i % 1000 == 9999) save();
	}
	
	std::cout << std::endl;
	// save();
	
	// save the image after tone mapping and gamma correction
	FILE* f;
	fopen_s(&f, "D:/learn/graphics/image.ppm", "w");
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i< w * h; i++) {
		fprintf(f, "%d %d %d ", toInt(image[i].x), toInt(image[i].y), toInt(image[i].z));
	}
	
}

void Scene::initScene() {
	prim = new Primitive*[20];
	numPrim = 9;
	prim[0] = new Plane(Vec(1, 0, 0), -1, Material(Color(0.75, 0.25, 0.25), 1.0, 0.0, 0.0));
	prim[1] = new Plane(Vec(-1, 0, 0), 99, Material(Color(0.25, 0.25, 0.25), 1.0, 0.0, 0.0));
	prim[2] = new Plane(Vec(0, 0, 1), 0, Material(Color(0.75, 0.75, 0.75), 1.0, 1.0, 0.0));
	prim[3] = new Plane(Vec(0, 0, -1), 170, Material(Color(0, 0, 0), 1.0, 0.0, 0.0));
	prim[4] = new Plane(Vec(0, 1, 0), 0, Material(Color(0.75, 0.75, 0.75), 1.0, 0.0, 0.0));
	prim[5] = new Plane(Vec(0, -1, 0), 81.6, Material(Color(0.75, 0.75, 0.75), 1.0, 0.0, 0.0));
	/*
	prim[0] = new Sphere(Vec(1e5 + 1, 40.8, 81.6), 1e5, Material(Color(.75, .25, .25), 1., 0., 0.));
	prim[1] = new Sphere(Vec(-1e5 + 99, 40.8, 81.6), 1e5, Material(Color(.25, .25, .25), 1., 0., 0.));
	prim[2] = new Sphere(Vec(50, 40.8, 1e5), 1e5, Material(Color(.75, .75, .75), 1., 0., 0.));
	prim[3] = new Sphere(Vec(50, 40.8, -1e5 + 170), 1e5, Material(Color(), 1., 0., 0.));
	prim[4] = new Sphere(Vec(50, 1e5, 81.6), 1e5, Material(Color(.75, .75, .75), 1., 0., 0.));
	prim[5] = new Sphere(Vec(50, -1e5 + 81.6, 81.6), 1e5, Material(Color(.75, .75, .75), 1., 0., 0.));
	*/

	prim[6] = new Sphere(Vec(27, 16.5, 47), 16.5, Material(Color(1, 1, 1)*0.999, 0.0, 1.0, 0.0));
	prim[7] = new Sphere(Vec(73, 16.5, 88), 16.5, Material(Color(1, 1, 1)*0.999, 0.0, 0.0, 1.0));
	prim[7]->material.SetRefrIndex(1.5);
	prim[8] = new Sphere(Vec(50, 8.5, 60), 8.5, Material(Color(1, 1, 1)*0.999, 1.0, 0.0, 0.0));
	// prim[0] = new Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(.75, .25, .25), DIFF),//Left
	// prim[1] = new Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(.25, .25, .75), DIFF),//Right
	// prim[2] = new Sphere(1e5, Vec(50, 40.8, 1e5), Vec(.75, .75, .75), DIFF),//Back
	// prim[3] = new Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), DIFF),//Front
	// prim[4] = new Sphere(1e5, Vec(50, 1e5, 81.6), Vec(.75, .75, .75), DIFF),//Bottomm
	// prim[5] = new Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(.75, .75, .75), DIFF),//Top
	// prim[6] = new Sphere(16.5, Vec(27, 16.5, 47), Vec(1, 1, 1)*.999, SPEC),//Mirror
	// prim[7] = new Sphere(16.5, Vec(73, 16.5, 88), Vec(1, 1, 1)*.999, REFR),//Glass
	// prim[8] = new Sphere(8.5, Vec(50, 8.5, 60), Vec(1, 1, 1)*.999, DIFF)
	
}
