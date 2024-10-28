#include "Cloth.h"

Cloth::Cloth(vec3 f_external)
{
	width = 40;
	height = 40;
	this->f_external = f_external;
	substep = 5;
	neighbors.resize((width - 1) * (height - 1) * 2);

	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			particle p;
			p.inverseMass = 0.0f;
			/*if (j == height - 1 && (i == 0 || i == width - 1))
				p.inverseMass = 0;*/
			p.p = vec3(0.2 * i, 0.2 * j, 0);
			p.x = p.p;
			p.v = vec3(0, 0, 0);
			particles.push_back(p);
		}
	}
	//
	for (int j = 0; j < height - 1; j++)
	{
		for (int i = 0; i < width - 1; i++)
		{
			vec3 p1, p2, p3;
			triangle t1,t2;
			t1.indices[0] = j * width + i;
			t1.indices[1] = j * width + i + 1;
			t1.indices[2] = j * width + i + width;
			t1.triangleIndex = triangles.size();
			triangles.push_back(t1);
			p1 = particles[t1.indices[0]].x;
			p2 = particles[t1.indices[1]].x;
			p3 = particles[t1.indices[2]].x;
			float mass = 1.0 /(0.5 * length(cross(p1 - p2, p1 - p3))) / 3;
			particles[t1.indices[0]].inverseMass += mass;
			particles[t1.indices[1]].inverseMass += mass;
			particles[t1.indices[2]].inverseMass += mass;



			t2.indices[0] = j * width + i + 1;
			t2.indices[1] = j * width + i + width;
			t2.indices[2] = j * width + i + width + 1;
			t2.triangleIndex = triangles.size();
			triangles.push_back(t2);
			p1 = particles[t2.indices[0]].x;
			p2 = particles[t2.indices[1]].x;
			p3 = particles[t2.indices[2]].x;
			mass = 1.0 / (0.5 * length(cross(p1 - p2, p1 - p3))) / 3;
			particles[t2.indices[0]].inverseMass += mass;
			particles[t2.indices[1]].inverseMass += mass;
			particles[t2.indices[2]].inverseMass += mass;

			addEdge(t1);
			addEdge(t2);
		}
	}
	particles[width * height - 1].inverseMass = 0.0;
	particles[width * (height-1)].inverseMass = 0.0;

}
float Cloth::calculatePhi(edge e, int triangleIndex)
{
	float phi;
	vec3 p1, p2, p3, p4;
	p1 = particles[e.indices[0]].x;
	p2 = particles[e.indices[1]].x;
	for (int i = 0; i < 3; i++)
	{
		int index_p3 = triangles[e.triangleIndex].indices[i];
		if (index_p3 != e.indices[0] && index_p3 != e.indices[1])
			p3 = particles[index_p3].x;
		int index_p4 = triangles[triangleIndex].indices[i];
		if (index_p4 != e.indices[0] && index_p4 != e.indices[1])
			p4 = particles[index_p4].x;
	}
	vec3 n1 = normalize(cross(p2 - p1, p3 - p1));
	vec3 n2 = normalize(cross(p2 - p1, p4 - p1));
	phi = acos(dot(n1, n2));
	return phi;
}
void Cloth::addEdge(triangle t)
{
	edge e1, e2, e3;
	e1.indices[0] = t.indices[0];
	e1.indices[1] = t.indices[1];
	e1.resLen = length(particles[e1.indices[0]].x - particles[e1.indices[1]].x);
	e1.triangleIndex = t.triangleIndex;

	e2.indices[0] = t.indices[0];
	e2.indices[1] = t.indices[2];
	e2.resLen = length(particles[e2.indices[0]].x - particles[e2.indices[1]].x);
	e2.triangleIndex = t.triangleIndex;

	e3.indices[0] = t.indices[1];
	e3.indices[1] = t.indices[2];
	e3.resLen = length(particles[e3.indices[0]].x - particles[e3.indices[1]].x);
	e3.triangleIndex = t.triangleIndex;

	if(edges.find(e1)==edges.end())
		edges.insert(e1);
	else
	{
		auto e = edges.find(e1);
		float phi = calculatePhi(e1,e->triangleIndex);
		neighbor n;
		n.phi = phi;
		n.triangleIndex = e1.triangleIndex;
		n.sharedEdge = e1;
		for (int i = 0; i < 3; i++)
		{
			int index_p3 = triangles[e1.triangleIndex].indices[i];
			if (index_p3 != e1.indices[0] && index_p3 != e1.indices[1])
				n.bendEdge[0] = index_p3;
			int index_p4 = triangles[e->triangleIndex].indices[i];
			if (index_p4 != e1.indices[0] && index_p4 != e1.indices[1])
				n.bendEdge[1] = index_p4;
		}
		n.resLen = length(particles[n.bendEdge[0]].x - particles[n.bendEdge[1]].x);
		//vec3 p1, p2, p3,p4;
		////wip
		//p1 = particles[e1.indices[0]].x;
		//p2 = particles[e1.indices[1]].x;
		//for (int i = 0; i < 3; i++)
		//{
		//	int index_p3 = triangles[e1.triangleIndex].indices[i];
		//	if(index_p3 != e1.indices[0] && index_p3 != e1.indices[1])
		//		p3 = particles[index_p3].x;
		//	int index_p4 = triangles[e->triangleIndex].indices[i];
		//	if (index_p4 != e1.indices[0] && index_p4 != e1.indices[1])
		//		p4 = particles[index_p4].x;
		//}
		//vec3 n1 = normalize(cross(p2-p1,p3-p1));
		//vec3 n2 = normalize(cross(p2-p1,p4-p1));
		//phi = acos(dot(n1,n2));
		neighbors[e->triangleIndex].push_back(n);
	}

	if (edges.find(e2) == edges.end())
		edges.insert(e2);
	else
	{
		auto e = edges.find(e2);
		float phi = calculatePhi(e2, e->triangleIndex);
		neighbor n;
		n.phi = phi;
		n.triangleIndex = e2.triangleIndex;
		n.sharedEdge = e2;
		for (int i = 0; i < 3; i++)
		{
			int index_p3 = triangles[e1.triangleIndex].indices[i];
			if (index_p3 != e2.indices[0] && index_p3 != e2.indices[1])
				n.bendEdge[0] = index_p3;
			int index_p4 = triangles[e->triangleIndex].indices[i];
			if (index_p4 != e2.indices[0] && index_p4 != e2.indices[1])
				n.bendEdge[1] = index_p4;
		}
		n.resLen = length(particles[n.bendEdge[0]].x - particles[n.bendEdge[1]].x);
		//vec3 p1, p2, p3, p4;
		////wip
		//p1 = particles[e2.indices[0]].x;
		//p2 = particles[e2.indices[1]].x;
		//for (int i = 0; i < 3; i++)
		//{
		//	int index_p3 = triangles[e2.triangleIndex].indices[i];
		//	if (index_p3 != e2.indices[0] && index_p3 != e2.indices[1])
		//		p3 = particles[index_p3].x;
		//	int index_p4 = triangles[e->triangleIndex].indices[i];
		//	if (index_p4 != e2.indices[0] && index_p4 != e2.indices[1])
		//		p4 = particles[index_p4].x;
		//}
		//vec3 n1 = normalize(cross(p2 - p1, p3 - p1));
		//vec3 n2 = normalize(cross(p2 - p1, p4 - p1));
		//phi = acos(dot(n1, n2));
		neighbors[e->triangleIndex].push_back(n);
	}

	if (edges.find(e3) == edges.end())
		edges.insert(e3);
	else
	{
		auto e = edges.find(e3);
		float phi = calculatePhi(e3, e->triangleIndex);
		neighbor n;
		n.phi = phi;
		n.triangleIndex = e3.triangleIndex;
		n.sharedEdge = e3;
		for (int i = 0; i < 3; i++)
		{
			int index_p3 = triangles[e1.triangleIndex].indices[i];
			if (index_p3 != e3.indices[0] && index_p3 != e3.indices[1])
				n.bendEdge[0] = index_p3;
			int index_p4 = triangles[e->triangleIndex].indices[i];
			if (index_p4 != e3.indices[0] && index_p4 != e3.indices[1])
				n.bendEdge[1] = index_p4;
		}
		n.resLen = length(particles[n.bendEdge[0]].x - particles[n.bendEdge[1]].x);
		//vec3 p1, p2, p3, p4;
		////wip
		//p1 = particles[e3.indices[0]].x;
		//p2 = particles[e3.indices[1]].x;
		//for (int i = 0; i < 3; i++)
		//{
		//	int index_p3 = triangles[e3.triangleIndex].indices[i];
		//	if (index_p3 != e3.indices[0] && index_p3 != e3.indices[1])
		//		p3 = particles[index_p3].x;
		//	int index_p4 = triangles[e->triangleIndex].indices[i];
		//	if (index_p4 != e3.indices[0] && index_p4 != e3.indices[1])
		//		p4 = particles[index_p4].x;
		//}
		//vec3 n1 = normalize(cross(p2 - p1, p3 - p1));
		//vec3 n2 = normalize(cross(p2 - p1, p4 - p1));
		//phi = acos(dot(n1, n2));
		neighbors[e->triangleIndex].push_back(n);
	}
}
void Cloth::preSolve(float dt)
{
	int index = -1;
	for (auto& p : particles)
	{
		index++;
		if (p.inverseMass == 0)
			continue;
		//p.v += dt * f_external;
		p.p = p.x;
		p.x += p.v * dt + p.inverseMass * dt * dt * f_external * 0.1f;
		
		//p.v *= 0.5;
	}
}


void Cloth::solve(float dt)
{
	solveEdges(dt);
	solveBending(dt);
}
void Cloth::solveEdges(float dt)
{
	for (auto& e : edges)
	{
		vec3 p1 = particles[e.indices[0]].x;
		vec3 p2 = particles[e.indices[1]].x;
		vec3 gradient = p2 - p1;
		float w1 = particles[e.indices[0]].inverseMass;
		float w2 = particles[e.indices[1]].inverseMass;
		float w = w1 + w2;
		float C = length(gradient) - e.resLen;
		float s = -C / (w + 0.001 / dt / dt);
		gradient = normalize(gradient);
		//if(!particles[e.indices[0]].pinned)
			particles[e.indices[0]].x -= s * w1 * gradient;
		//if (!particles[e.indices[1]].pinned)
			particles[e.indices[1]].x += s * w2 * gradient;

	}
}
void Cloth::solveBending(float dt)
{
	for (int i = 0; i < neighbors.size(); i++)
	{
		for (int j = 0; j < neighbors[i].size(); j++)
		{
			/*float beta = calculatePhi(neighbors[i][j].sharedEdge, i);
			float C = beta - neighbors[i][j].phi;
			float s = i;*/
			vec3 p1 = particles[neighbors[i][j].bendEdge[0]].x;
			vec3 p2 = particles[neighbors[i][j].bendEdge[1]].x;
			vec3 gradient = p2 - p1;
			float w1 = particles[neighbors[i][j].bendEdge[0]].inverseMass;
			float w2 = particles[neighbors[i][j].bendEdge[1]].inverseMass;
			float w = w1 + w2;
			float C = length(gradient) - neighbors[i][j].resLen;
			float s = -C / (w + 0.01 /dt/dt);
			gradient = normalize(gradient);
			particles[neighbors[i][j].bendEdge[0]].x -= s * w1 * gradient;
			particles[neighbors[i][j].bendEdge[1]].x += s * w2 * gradient;

		}
	}



}
void Cloth::postSolve(float dt)
{
	for (auto& p : particles)
	{
		//if(!p.pinned)
			p.v = (p.x - p.p) / dt;
			if (length(p.v) > 10.0f)
				p.v = normalize(p.v) * 10.0f;
			//p.v *= 0.99;
	}
}


void Cloth::update(float dt)
{
	float sub_dt = dt / substep;
	for (int i = 0; i < substep; i++)
	{
		preSolve(sub_dt);
		solve(sub_dt);
		postSolve(sub_dt);
	}
}
void Cloth::draw()
{
	int n = particles.size() * 3;
	float* vertices = new float[n];
	for (int i = 0; i < particles.size(); i++)
	{
		vertices[3 * i] = particles[i].x.x;
		vertices[3 * i + 1] = particles[i].x.y;
		vertices[3 * i + 2] = particles[i].x.z;
	}
	int* indices = new int[triangles.size() * 3];
	for (int i = 0; i < triangles.size(); i++)
	{
		indices[3 * i] = triangles[i].indices[0];
		indices[3 * i + 1] = triangles[i].indices[1];
		indices[3 * i + 2] = triangles[i].indices[2];
	}

	GLuint VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * n, vertices, GL_STATIC_DRAW);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * triangles.size() * 3, indices, GL_STATIC_DRAW);

	
	glDrawElements(GL_TRIANGLES, triangles.size() * 3, GL_UNSIGNED_INT, 0);

	glBindVertexArray(0);

}