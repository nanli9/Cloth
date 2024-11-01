#pragma once
#ifndef SPATIALHASH_H
#define SPATIALHASH_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace glm;

class spatialHash
{
public:
	float spacing;
	int** table;

	spatialHash(int numOfObjects);
	int hashCoords(vec3 pos);
	int* query(vec3 pos, float d, int& querySize);
	vec3 intCoord(vec3 pos);

};


#endif // !SPATIALHASH_H
