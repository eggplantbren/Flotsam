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
	latent_boost = randomU();
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
		latent_boost += pow(10., 1.5 - 6.*randomU())*randn();
		latent_boost = mod(latent_boost, 1.);
	}

	return logH;
}

double NoiseProperties::get_boost() const
{
	if(latent_boost < 0.5)
		return 1.;

	double b = 2.*(latent_boost - 0.5);
	b = exp(log(1.) + log(100.)*b);
	return b;
}

void NoiseProperties::print(ostream& out) const
{
	out<<nu<<' '<<get_boost()<<' ';
}

