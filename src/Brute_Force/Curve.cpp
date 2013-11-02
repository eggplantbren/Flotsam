#include "Curve.h"
#include "RandomNumberGenerator.h"
#include <cassert>
#include <cmath>

using namespace std;
using namespace DNest3;

Curve::Curve(double t_min, double t_max, int N)
:N(N)
,t_min(t_min)
,t_max(t_max)
,t_range(t_max - t_min)
,n(N)
,y(N)
{
	assert(t_max > t_min);
}

void Curve::fromPrior()
{
	for(int i=0; i<N; i++)
		n[i] = randn();
	L = exp(log(1E-2*t_range) + log(1E4)*randomU());
	assemble();
}

double Curve::perturb()
{

	assemble();
}

void Curve::assemble()
{
	double a = exp(-1./L);
	double b = sqrt(1. - a*a);

	y[0] = n[0];
	for(int i=1; i<N; i++)
		y[i] = a*y[i-1] + b*n[i];
}

void Curve::print(ostream& out) const
{
	for(int i=0; i<N; i++)
		out<<y[i]<<' ';
}

