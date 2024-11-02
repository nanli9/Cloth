#include "Cloth.h"

Cloth::Cloth(vec3 f_external)
{
	width = 40;
	height = 50;
	this->f_external = f_external;
	substep = 5;
	neighbors.resize((width - 1) * (height - 1) * 2);
	stretchCompliance = 0.001f;
	bendCompliance = 0.001f;
	grabIndex = -1;
	totalMass = 0;
	k_damping = 0.0f;
	lineDisplay = false;
	pinned = false;
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			particle p;
			p.inverseMass = 0.0f;
			p.p = vec3(0.1 * i, 0 ,0.1 * j);
			p.restPos = p.p;
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
			float mass = 1.0 / (0.5 * length(cross(p1 - p2, p1 - p3))) / 3;
			particles[t1.indices[0]].inverseMass += mass;
			particles[t1.indices[1]].inverseMass += mass;
			particles[t1.indices[2]].inverseMass += mass;



			t2.indices[0] = j * width + i + 1;
			t2.indices[1] = j * width + i + width + 1;
			t2.indices[2] = j * width + i + width;
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
	for (auto& p : particles)
		totalMass += 1.0f / p.inverseMass;
	

	this->hash = new spatialHash(particles.size(), 0.1f);
	
	pinnedParticleInverseMass[0] = particles[width * height - 1].inverseMass;
	pinnedParticleInverseMass[1] = particles[width * (height - 1)].inverseMass;
	pinnedParticleInverseMass[2] = particles[0].inverseMass;
	pinnedParticleInverseMass[3] = particles[width - 1].inverseMass;
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
void Cloth::preSolve(float dt, vector<RigidBody>& rigidBodies)
{
	vec3 xcm(0,0,0), vcm(0, 0, 0),L(0,0,0);
	mat3 I(0.0f);
	for (auto& p : particles)
	{
		if (p.inverseMass == 0)
			continue;
		xcm += p.x / p.inverseMass;
		vcm += p.v / p.inverseMass;
	}
	xcm /= totalMass;
	vcm /= totalMass;
	for (auto& p : particles)
	{
		if (p.inverseMass == 0)
			continue;
		vec3 r = p.x - xcm;
		L += cross(r, p.v / p.inverseMass);
		mat3 r_star(0.0f);
		r_star[1][0] = -r.z;
		r_star[2][0] = r.y;
		r_star[0][1] = r.z;
		r_star[2][1] = -r.x;
		r_star[0][2] = -r.y;
		r_star[1][2] = r.x;
		I += r_star * transpose(r_star) * (1.0f / p.inverseMass);
	}
	vec3 omega = inverse(I) * L;
	for (auto& p : particles)
	{
		if (p.inverseMass == 0)
			continue;
		p.p = p.x;
		p.v += p.inverseMass * dt * f_external * 0.1f;
		//damping
		vec3 delta_v = vcm + cross(omega, p.x - xcm) - p.v;
		p.v += k_damping * delta_v;
		p.x += p.v * dt;
		if (p.x.y <= -2.0f + 0.5 * 0.01)
		{
			vec3 v = (p.x - p.p);
			p.x -= v;
			p.x.y = -2.0f + 0.5 * 0.01;

		}
		handleCollision(p.x, rigidBodies);
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
		float s = -C / (w + stretchCompliance / dt / dt);
		gradient = normalize(gradient);
		particles[e.indices[0]].x -= s * w1 * gradient;
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
			float s = -C / (w + bendCompliance /dt/dt);
			gradient = normalize(gradient);
			particles[neighbors[i][j].bendEdge[0]].x -= s * w1 * gradient;
			particles[neighbors[i][j].bendEdge[1]].x += s * w2 * gradient;

		}
	}

}
void Cloth::selfCollision(float dt)
{
	for (int i = 0; i < particles.size(); i++)
	{
		float d = 0.0f;
		hash->query(particles[i].x, 0.01 ,i);
		int x = hash->querySize;
		for (int j = 0; j < hash->querySize; j++)
		{
			int index_1 = i;
			int index_2 = hash->queryIndices[j];

			vec3 p = particles[index_1].x - particles[index_2].x;
			vec3 p_0 = particles[index_1].restPos - particles[index_2].restPos;
			float restLen2 = dot(p_0, p_0);
			if (dot(p, p) <= std::min(0.0001f, restLen2))
			{
				/*float l = sqrt(dot(p, p));
				particles[index_1].x -= 0.5f * p * (0.01f - l) / l;
				particles[index_2].x += 0.5f * p * (0.01f - l) / l;*/

				vec3 v1 = (particles[index_1].x - particles[index_1].p) / dt;
				vec3 v2 = (particles[index_2].x - particles[index_2].p) / dt;
				vec3 v_avg = (v1 + v2) / 2.0f;
				particles[index_1].v = v1;
				particles[index_2].v = v2;
				particles[index_1].x += d * (v_avg - v1) * dt;
				particles[index_2].x += d * (v_avg - v2) * dt;
			}
		}
	}
}
void Cloth::postSolve(float dt)
{
	for (auto& p : particles)
	{
		p.v = (p.x - p.p) / dt;
		if (length(p.v) > 10.0f)
			p.v = normalize(p.v) * 10.0f;
	}
}


void Cloth::update(float dt, vector<RigidBody>& rigidBodies)
{
	if (pinned)
	{
		particles[width * height - 1].inverseMass = 0.0;
		particles[width * (height-1)].inverseMass = 0.0;
		particles[0].inverseMass = 0.0;
		particles[width-1].inverseMass = 0.0;
	}
	hash->clean();
	for (int i = 0;i<particles.size();i++)
	{
		//hash->table.insert(pair<vec3, int>(particles[i].x,i));
		hash->insert(particles[i].x, i);
	}
	float sub_dt = dt / substep;
	for (int i = 0; i < substep; i++)
	{
		preSolve(sub_dt, rigidBodies);
		solve(sub_dt);
		selfCollision(sub_dt);
		postSolve(sub_dt);
	}
}
void Cloth::draw(const Shader& shader, mat4& shadowMatrix)
{
	int n = particles.size() * 6;
	float* vertices = new float[n];
	for (int i = 0; i < particles.size(); i++)
	{
		vertices[6 * i] = particles[i].x.x;
		vertices[6 * i + 1] = particles[i].x.y;
		vertices[6 * i + 2] = particles[i].x.z;
		vertices[6 * i + 3] = 0.0f;
		vertices[6 * i + 4] = 0.0f;
		vertices[6 * i + 5] = 0.0f;
	}
	int* indices = new int[triangles.size() * 3];
	for (int i = 0; i < triangles.size(); i++)
	{
		indices[3 * i] = triangles[i].indices[0];
		indices[3 * i + 1] = triangles[i].indices[1];
		indices[3 * i + 2] = triangles[i].indices[2];
		vec3 v1 = particles[indices[3 * i]].x - particles[indices[3 * i + 1]].x;
		vec3 v2 = particles[indices[3 * i]].x - particles[indices[3 * i + 2]].x;
		vec3 n = normalize(cross(v2, v1));

		vertices[6 * indices[3 * i] + 3] += n.x;
		vertices[6 * indices[3 * i] + 4] += n.y;
		vertices[6 * indices[3 * i] + 5] += n.z;

		vertices[6 * indices[3 * i + 1] + 3] += n.x;
		vertices[6 * indices[3 * i + 1] + 4] += n.y;
		vertices[6 * indices[3 * i + 1] + 5] += n.z;

		vertices[6 * indices[3 * i + 2] + 3] += n.x;
		vertices[6 * indices[3 * i + 2] + 4] += n.y;
		vertices[6 * indices[3 * i + 2] + 5] += n.z;
	}

	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);
	
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(float) * n, vertices, GL_DYNAMIC_DRAW);

	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * triangles.size() * 3, indices, GL_DYNAMIC_DRAW);

	if (lineDisplay)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(2);
	}

	glDrawElements(GL_TRIANGLES, triangles.size() * 3, GL_UNSIGNED_INT, 0);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	shader.setMat4("model", shadowMatrix);
	shader.setInt("shadow", 1);
	glDrawElements(GL_TRIANGLES, triangles.size() * 3, GL_UNSIGNED_INT, 0);
	shader.setInt("shadow", 0);

	glBindVertexArray(0);

}
void Cloth::grab(vec3 pos)
{
	float d = 0.1;
	int index = -1;
	for (int i = 0; i < particles.size(); i++)
	{
		float d_new = length(particles[i].x - pos);
		if (d_new < d)
		{
			d = d_new;
			index = i;
		}

	}
	grabIndex = index;
	if (grabIndex != -1)
	{
		grabPointInvereMass = particles[grabIndex].inverseMass;
		particles[grabIndex].inverseMass = 0.0f;
	}
}
void Cloth::handleCollision(vec3& p, vector<RigidBody>& rigidBodies)
{
	vec3 center = vec3(1.5,-1,2);
	float r = 1;
	vec3 v1 = p - center;
	if (dot(v1, v1) <= r * r)
	{
		v1 = r * normalize(v1);
		p = v1 + center;
	}

}
void Cloth::reset()
{
	pinned = false;
	particles[width * height - 1].inverseMass = pinnedParticleInverseMass[0];
	particles[width * (height - 1)].inverseMass = pinnedParticleInverseMass[1];
	particles[0].inverseMass = pinnedParticleInverseMass[2];
	particles[width - 1].inverseMass = pinnedParticleInverseMass[3];
	f_external = vec3(0, 0, 0);
	for (int j = 0; j < height; j++)
	{
		for (int i = 0; i < width; i++)
		{
			particles[j * width + i].x = particles[j * width + i].restPos;
			particles[j * width + i].p = particles[j * width + i].restPos;
			particles[j * width + i].v = vec3(0,0,0);
		}
	}
}