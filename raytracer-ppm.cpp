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


int main(int argc, char *argv[]) {
	Engine eng;
	eng.render();
	system("pause");
}