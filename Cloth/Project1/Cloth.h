#pragma once
#ifndef CLOTH_H
#define CLOTH_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <set>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace glm;
using namespace std;

struct particle
{
	vec3 v; // velocity
	vec3 p; //previos pos
	vec3 x; //cur pos
	float inverseMass;
	//bool pinned;
};
struct edge
{
	int indices[2];
	float resLen;
	int triangleIndex;
	bool operator<(const edge& other) const {
		if(indices[0] == other.indices[0])
			return indices[1] < other.indices[1];
		return indices[0] < other.indices[0];
	}
};
struct neighbor
{
	float phi;
	int triangleIndex;
	edge sharedEdge;
	int PointsIndices[4];
	int bendEdge[2];
	float resLen;
};
struct triangle
{
	int indices[3];
	int triangleIndex;
};
class Cloth
{
public:
	vector<particle> particles;
	vector<triangle> triangles;
	vector<vector<neighbor>> neighbors;
	float bendCompliance;
	float stretchCompliance;
	set<edge> edges;
	int height, width;
	vec3 f_external;
	int substep;
	int grabIndex;
	int grabPointInvereMass;
	Cloth(vec3 f_external);
	void preSolve(float dt);
	void solve(float dt);
	void solveEdges(float dt);
	void solveBending(float dt);
	void postSolve(float dt);
	void update(float dt);
	void draw();
	void addEdge(triangle t);
	float calculatePhi(edge e,int triangleIndex);
	void grab(vec3 pos);
};


#endif // !CLOTH_H
