#include <cmath>
#include "NoiseProperties.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"

using namespace DNest3;
using namespace std;

NoiseProperties::NoiseProperties(int N)
:u(N)
{

}

void NoiseProperties::fromPrior()
{
	f_bad = randomU();
	boost1 = randomU();
	boost2 = randomU();

	for(size_t i=0; i<u.size(); i++)
		u[i] = randomU();
}

double NoiseProperties::perturb()
{
	double logH = 0.;

	int which = randInt(3);

	if(which == 0)
	{
		f_bad += pow(10., 1.5 - 6.*randomU())*randn();
		f_bad = mod(f_bad, 1.);
	}
	else if(which == 1)
	{
		boost1 += pow(10., 1.5 - 6.*randomU())*randn();
		boost1 = mod(boost1, 1.);
	}
	else if(which == 2)
	{
		boost2 += pow(10., 1.5 - 6.*randomU())*randn();
		boost2 = mod(boost2, 1.);
	}

	double chance = pow(10., 0.5 - 6.*randomU());
	for(size_t i=0; i<u.size(); i++)
		if(randomU() <= chance)
			u[i] = randomU();

	return logH;
}

double NoiseProperties::get_boost(int i) const
{
	if(i < 0 || i >= static_cast<int>(u.size()))
		return 1.;

	double b1 = (boost1 < 0.5)?(0.):((boost1 - 0.5)/0.5);
	double b2 = boost2;

	b1 = exp(log(1.) + log(100.)*b1);
	b2 = exp(log(1.) + log(100.)*b2);

	return (u[i] < f_bad)?(b1*b2):(b1);
}

