#include "spatialHash.h"


spatialHash::spatialHash(int numOfObjects, float spacing)
{
	this->tableSize = 2 * numOfObjects;
	queryIndices.resize(this->tableSize);
	querySize = 0;
	table.resize(this->tableSize);
	this->spacing = spacing;
	for (auto& i : table)
		i.reserve(4);
}
void spatialHash::insert(const vec3 p,int index)
{
	long long x = p.x / spacing;
	long long y = p.y / spacing;
	long long z = p.z / spacing;
	if (/*index == 0 ||*/ index == 477 || index == 1280)
	{
		int asd = index;
	}
	int h = hashCoords(x,y,z);
	table[h].push_back(index);
}


int spatialHash::hashCoords(long long x, long long y, long long z)
{
	long long hash = (x * 92837111) ^ (y * 689287499) ^ (z * 283923481);
	long long a = x * 92837111;
	return abs(hash) % tableSize;
}
void spatialHash::query(vec3 pos, float d, int particleIndex)
{
	long long x0 = (pos.x - d) / spacing;
	long long y0 = (pos.y - d) / spacing;
	long long z0 = (pos.z - d) / spacing;

	long long x1 = (pos.x + d) / spacing;
	long long y1 = (pos.y + d) / spacing;
	long long z1 = (pos.z + d) / spacing;

	querySize = 0;

	for (long long i = x0; i <= x1; i++)
	{
		for (long long j = y0; j <= y1; j++)
		{
			for (long long k = z0; k <= z1; k++)
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
