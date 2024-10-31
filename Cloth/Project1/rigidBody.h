#pragma once
#ifndef RIGID_BODY_H
#define RIGID_BODY_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <set>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace glm;
using namespace std;
struct Y
{
	vec3 x;
	mat3 R;
	vec3 P;
	vec3 L;

	mat3 I_0;
	mat3 I;
	mat3 I_inverse;

	vec3 omega;

	vec3 force;
	vec3 torque;
};

class RigidBody
{
public:
	Y y_t;
	float mass;
	unsigned int VBO, cubeVAO;
	vector<vec3> vertices;
	RigidBody();
	RigidBody(vec3 x, mat3 R, vec3 P, vec3 L,float mass);
	void draw();


};



#endif // !RIGID_BODY_H

