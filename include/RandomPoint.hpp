
#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class RandomPoint {
public:
	RandomPoint(int x, float y) :seed(x), random(y) {};
	RandomPoint() :seed(0), random(0.0f) {};

	void calcRandom(int& _seed) {
		const int m = pow(2, 31) - 1;
		const int a = 16807;
		const int q = 127773;
		const int r = 2836;

		int temp = a * (_seed % q) - r * (_seed / q);

		if (temp < 0)
		{
			temp = m + temp;
		}
		//num is the seed number
		_seed = temp;
		float t = _seed * 1.0f / m;

		this->random = t;
		this->seed = _seed;
	}
public:
	int seed;
	float random;
};

