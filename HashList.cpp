#include "stdafx.h"
#include "HashList.h"
#include <iostream>

void HashList::build_hash_grid(int w, int h){
	// find the bounding box of all the measurement points
	hpBox.reset();//包围了所有hitPoint的包围盒
	
	for (int i = 0; i < hList.size(); ++ i){
		hpBox.fit(hList[i]->pos);
	}
	
	// heuristic for initial radius
	Vec ssize = hpBox.max - hpBox.min;//包围盒的大小
	double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;//半径大小？？
	
	// determine hash table size
	// we now find the bounding box of all the measurement points inflated by the initial radius
	hpBox.reset();
	int vphoton = int(hList.size());
	for (int i = 0; i < hList.size(); ++ i){
		hList[i]->radius2 = irad * irad;
		hList[i]->N = 0;
		hList[i]->flux = Vec();
		hpBox.fit(hList[i]->pos-irad);
		hpBox.fit(hList[i]->pos+irad);
		
	}
	
	// make each grid cell two times larger than the initial radius
	hash_s = 1.0/(irad*2.0);
	num_hash = vphoton;
	
	// build the hash table
	hashGrid = new std::vector<HPoint*>[num_hash];
	// for (unsigned int i=0; i<num_hash;i++) hash_grid[i] = NULL;
	for(int i = 0; i < hList.size(); ++ i){
		Vec BMin = ((hList[i]->pos - irad) - hpBox.min) * hash_s;
		Vec BMax = ((hList[i]->pos + irad) - hpBox.min) * hash_s;
		
		for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++){
			for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++){
				for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++){
					int hv = hash(ix, iy, iz);
					hashGrid[hv].push_back(hList[i]);
				}
			}
		}
	}
	/*
	for (int i = 0; i < num_hash; ++i) {
		for (int j = 0; j < hashGrid[i].size(); ++j) {
			std::cout << "(" << hashGrid[i][j]->x << "," << hashGrid[i][j]->y << ")";
		}
		std::cout << std::endl;
	}
	*/
}
