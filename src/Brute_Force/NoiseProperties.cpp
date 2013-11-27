#include <cmath>
#include "NoiseProperties.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"

using namespace DNest3;
using namespace std;

NoiseProperties::NoiseProperties()
{

}

void NoiseProperties::fromPrior()
{
	nu = exp(log(0.5) + log(200.)*randomU());
	boost = exp(log(1.) + log(100.)*randomU());
}

double NoiseProperties::perturb()
{
	double logH = 0.;

	int which = randInt(2);

	if(which == 0)
	{
		nu = log(nu);
		nu += log(200.)*pow(10., 1.5 - 6.*randomU())*randn();
		nu = mod(nu - log(0.5), log(200.)) + log(0.5);
		nu = exp(nu);
	}
	else if(which == 1)
	{
		boost = log(boost);
		boost += log(100.)*pow(10., 1.5 - 6.*randomU())*randn();
		boost = mod(boost - log(1.), log(100.)) + log(1.);
		boost = exp(boost);
	}

	return logH;
}

void NoiseProperties::print(ostream& out) const
{
	out<<nu<<' '<<boost<<' ';
}

