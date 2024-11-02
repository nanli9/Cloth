#include "spatialHash.h"


spatialHash::spatialHash(int numOfObjects, float spacing)
{
	this->tableSize = numOfObjects;
	queryIndices.resize(this->tableSize);
	querySize = 0;
	table.resize(this->tableSize);
	this->spacing = spacing;
	for (auto& i : table)
		i.reserve(4);
}
void spatialHash::insert(const vec3 p,int index)
{
	int x = p.x / spacing;
	int y = p.y / spacing;
	int z = p.z / spacing;

	int h = hashCoords(x,y,z);
	table[h].push_back(index);
}


int spatialHash::hashCoords(int x, int y, int z)
{
	int hash = (x * 92837111) ^ (y * 689287499) ^ (z * 283923481);
	return abs(hash) % tableSize;
}
void spatialHash::query(vec3 pos, float d, int particleIndex)
{
	int x0 = (pos.x - d) / spacing;
	int y0 = (pos.y - d) / spacing;
	int z0 = (pos.z - d) / spacing;

	int x1 = (pos.x + d) / spacing;
	int y1 = (pos.y + d) / spacing;
	int z1 = (pos.z + d) / spacing;

	querySize = 0;

	for (int i = x0; i <= x1; i++)
	{
		for (int j = y0; j <= y1; j++)
		{
			for (int k = z0; k <= z1; k++)
			{
				int h = hashCoords(i,j,k);
				for (int x : table[h])
				{
					if (x != particleIndex)
					{
						queryIndices[querySize] = x;
						querySize++;
					}
				}
			}
		}
	}
}

void spatialHash::clean()
{
	querySize = 0;
	table.resize(this->tableSize);
	for (auto& i : table)
	{
		i.clear();
		i.reserve(4);
	}
}
