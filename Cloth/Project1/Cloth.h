#pragma once
#ifndef CLOTH_H
#define CLOTH_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <set>
#include <iostream>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "rigidBody.h"
#include "spatialHash.h"
#include "Shader.h"

using namespace glm;
using namespace std;

struct particle
{
	vec3 v; // velocity
	vec3 p; //previos pos
	vec3 x; //cur pos
	float inverseMass;
	vec3 restPos;
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
	spatialHash* hash;
	vector<particle> particles;
	vector<triangle> triangles;
	vector<vector<neighbor>> neighbors;
	float pinnedParticleInverseMass[4];
	float bendCompliance;
	float stretchCompliance;
	bool lineDisplay;
	set<edge> edges;
	int height, width;
	vec3 f_external;
	int substep;
	int grabIndex;
	int grabPointInvereMass;
	float totalMass;
	float k_damping;
	bool pinned;
	Cloth(vec3 f_external);
	void preSolve(float dt, vector<RigidBody>& rigidBodies);
	void solve(float dt);
	void solveEdges(float dt);
	void solveBending(float dt);
	void postSolve(float dt);
	void update(float dt, vector<RigidBody>& rigidBodies);
	void draw(const Shader& shader, mat4& shadowMatrix);
	void addEdge(triangle t);
	float calculatePhi(edge e,int triangleIndex);
	void grab(vec3 pos);
	void handleCollision(vec3& p,vector<RigidBody>& rigidBodies);
	void selfCollision(float dt);
	void reset();
private:
	GLuint VAO, VBO, EBO;
};


#endif // !CLOTH_H
