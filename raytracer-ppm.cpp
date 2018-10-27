// raytracer-ppm.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
// expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)

#include "Vec.h"
#include "Primitive.h"

#include <math.h>   // smallppm, Progressive Photon Mapping by T. Hachisuka
#include <stdlib.h> // originally smallpt, a path tracer by Kevin Beason, 2008
#include <stdio.h>  // Usage: ./smallppm 100000 && xv image.ppm
#include "Engine.h"

double anotherToDouble(double x) {
	return pow(1 - exp(-x), 1 / 2.2);
}

int main(int argc, char *argv[]) {
	Engine* eng;
	int num = 50;
	eng = new Engine[num];
	srand(19981102);
	
	int w = 1024;
	int h = 768;

	#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < num; ++i) {
		std::cout << i << std::endl;
		// Engine eng;
		eng[i].render(i);
	}
	Bmp img(h, w);
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w; ++j) {
			Color tmp(0, 0, 0);
			for (int b = 0; b < num; ++b) {
				tmp = tmp + eng[b].image[i * w + j];
			}
			tmp = tmp * (1. / num);
			// if (tmp.x > 1) tmp.x = 1;
			// if (tmp.y > 1) tmp.y = 1;
			// if (tmp.z > 1) tmp.z = 1;
			img.SetColor(i, j, Color(anotherToDouble(tmp.x), anotherToDouble(tmp.y), anotherToDouble(tmp.z)));
			// bp.SetColor(i, j, Color(toDouble(image[i * w + j].x), toDouble(image[i * w + j].y), toDouble(image[i * w + j].z)));
		}
	}
	img.Output("D:/learn/graphics/image_sppm.bmp");
	

	system("pause");
}