#pragma once
#ifndef SPATIALHASH_H
#define SPATIALHASH_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <unordered_map>

using namespace glm;
using namespace std;


class spatialHash
{
public:
	float spacing;
	int tableSize;
	vector<int> queryIndices;
	vector<vector<int>> table;
	int querySize;
	spatialHash(int numOfObjects,float spacing);
	int hashCoords(int x, int y, int z);
	void query(vec3 pos, float d,int particleIndex);
	void insert(const vec3 p, int index);
	void clean();
	//vec3 intCoord(vec3 pos);
};


#endif // !SPATIALHASH_H
