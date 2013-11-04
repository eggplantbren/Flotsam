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

	mu = 10.*tan(M_PI*(randomU() - 0.5));
	beta = exp(log(1E-3) + log(1E6)*randomU());
	L = exp(log(1E-2*t_range) + log(1E4)*randomU());
	assemble();
}

double Curve::perturb()
{
	double logH = 0.;
	if(randomU() <= 0.5)
		logH += perturb1();
	else
		logH += perturb2();
	return logH;
}

double Curve::log_prob() const
{
	double logP = 0.;

	double a = exp(-1./L);
	double sig0 = beta/sqrt(1. - a*a);
	logP += -log(sig0) - 0.5*pow((y[0] - mu)/sig0, 2);
	for(int i=1; i<N; i++)
		logP += -log(beta) - 0.5*pow((y[i] - a*y[i-1])/beta, 2);

	return logP;
}

double Curve::perturb_param()
{
	int which = randInt(3);

	if(which == 0)
	{
		mu = atan(mu/10.)/M_PI;
		mu += pow(10., 1.5 - 6.*randomU())*randn();
		mu = mod(mu, 1.);
		mu = 10.*tan(M_PI*(mu - 0.5));
	}
	else if(which == 1)
	{
		beta = log(beta);
		beta += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		beta = mod(beta - log(1E-3), log(1E6)) + log(1E-3);
		beta = exp(beta);
	}
	else if(which == 2)
	{
		L = log(L);
		L += log(1E4)*pow(10., 1.5 - 6.*randomU())*randn();
		L = mod(L - log(1E-2*t_range), log(1E4)) + log(1E-2*t_range);
		L = exp(L);
	}
	return 0.;
}

double Curve::perturb2()
{
	double logH = 0.;

	logH -= log_prob();
	logH += perturb_param();
	logH += log_prob();
	disassemble();

	return logH;
}

double Curve::perturb1()
{
	double logH = 0.;

	logH += perturb_param();

	double chance = pow(10., 0.5 - 6.*randomU());
	for(int i=0; i<N; i++)
		if(randomU() <= chance)
			n[i] = randn();

	assemble();
	return logH;
}

void Curve::assemble()
{
	double a = exp(-1./L);
	y[0] = n[0]*beta/sqrt(1. - a*a);
	for(int i=1; i<N; i++)
		y[i] = a*y[i-1] + beta*n[i];
}

void Curve::disassemble()
{
	double a = exp(-1./L);
	n[0] = y[0]*sqrt(1. - a*a)/beta;
	for(int i=1; i<N; i++)
		n[i] = (y[i] - a*y[i-1])/beta;
}


double Curve::evaluate(double t) const
{
	int i = static_cast<int>((t - t_min)/dt);
	double w = (t - (t_min + i*dt))/dt;
	if(i >= 0 && i < N - 1)
		return mu + (1. - w)*y[i] + w*y[i+1];
	return mu;
}

void Curve::print(ostream& out) const
{
	for(int i=0; i<N; i++)
		out<<y[i]<<' ';
}

