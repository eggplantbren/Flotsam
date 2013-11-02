#include "Curve.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cassert>
#include <cmath>

using namespace std;
using namespace DNest3;

Curve::Curve(double t_min, double t_max, int N)
:N(N)
,t_min(t_min)
,t_max(t_max)
,t_range(t_max - t_min)
,dt(t_range/(N - 1))
,n(N)
,y(N)
{
	assert(t_max > t_min);
}

void Curve::fromPrior()
{
	for(int i=0; i<N; i++)
		n[i] = randn();

	beta = exp(log(1E-3) + log(1E6)*randomU());
	L = exp(log(1E-2*t_range) + log(1E4)*randomU());
	assemble();
}

double Curve::perturb()
{
	int which = randInt(3);

	if(which == 0)
	{
		beta = log(beta);
		beta += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		beta = mod(beta - log(1E-3), log(1E6)) + log(1E-3);
		beta = exp(beta);
	}
	else if(which == 1)
	{
		L = log(L);
		L += log(1E4)*pow(10., 1.5 - 6.*randomU())*randn();
		L = mod(L - log(1E-2*t_range), log(1E4)) + log(1E-2*t_range);
		L = exp(L);
	}
	else
	{
		double chance = pow(10., 0.5 - 6.*randomU());
		for(int i=0; i<N; i++)
			if(randomU() <= chance)
				n[i] = randn();
	}

	assemble();
	return 0.;
}

void Curve::assemble()
{
	double a = exp(-1./L);

	y[0] = n[0]*beta/sqrt(1. - a*a);
	for(int i=1; i<N; i++)
		y[i] = a*y[i-1] + beta*n[i];
}

double Curve::evaluate(double t) const
{
	int i = static_cast<int>((t - t_min)/dt);
	double w = (t - (t_min + i*dt))/dt;
	if(i >= 0 && i < N - 1)
		return (1. - w)*y[i] + w*y[i+1];
	return 0.;
}

void Curve::print(ostream& out) const
{
	for(int i=0; i<N; i++)
		out<<y[i]<<' ';
}

