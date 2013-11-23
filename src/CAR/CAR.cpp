#include "CAR.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include "RandomNumberGenerator.h"
#include "Utils.h"

using namespace std;
using namespace DNest3;

CAR::CAR(const vector<double>& times)
:t(times)
,n(times.size())
,y(times.size())
{
	double t_range = *max_element(times.begin(), times.end()) -
				*min_element(times.begin(), times.end());
	L_min = 1E-2*t_range;
	L_max = 1E2 *t_range;
}

void CAR::fromPrior()
{
	for(size_t i=0; i<n.size(); i++)
		n[i] = randn();
		
	mu = 10.*tan(M_PI*(randomU() - 0.5));
	sigma = exp(log(1E-3) + log(1E6)*randomU());
	L = exp(log(L_min) + log(L_max/L_min)*randomU());

	assemble();
}

double CAR::perturb()
{
	double logH = 0.;

	int which = randInt(3);
	if(which == 0)
	{
		mu = 0.5 + atan(mu/10.)/M_PI;
		mu += pow(10., 1.5 - 6.*randomU())*randn();
		mu = mod(mu, 1.);
		mu = 10.*tan(M_PI*(mu - 0.5));
	}
	else if(which == 1)
	{
		sigma = log(sigma);
		sigma += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		sigma = mod(sigma - log(1E-3), log(1E6)) + log(1E-3);
		sigma = exp(sigma);
	}
	else if(which == 2)
	{
		L = log(L);
		L += log(L_max/L_min)*pow(10., 1.5 - 6.*randomU())*randn();
		L = mod(L - log(L_min), log(L_max/L_min)) + log(L_min);
		L = exp(L);
	}

	// Now do the ns
	if(randomU() <= 0.5)
	{
		which = randInt(n.size());
		logH -= -0.5*pow(n[which], 2);
		n[which] += pow(10., 1.5 - 6.*randomU())*randn();
		logH += -0.5*pow(n[which], 2);
	}
	else if(which == 1)
	{
		double chance = pow(10., 0.5 - 4.*randomU());
		for(size_t i=0; i<n.size(); i++)
			if(randomU() <= chance)
				n[i] = randn();
	}

	assemble();
	return logH;
}

void CAR::assemble()
{
	y[0] = mu + sigma*n[0];

	double gap, mean, sd;
	for(size_t i=1; i<n.size(); i++)
	{
		gap = t[i] - t[i-1];
		mean = mu + (y[i-1] - mu)*exp(-gap/L);
		sd = sigma*sqrt(1. - exp(-2.*gap/L));
		y[i] = mean + sd*n[i];
	}
}

void CAR::print(ostream& out) const
{
	out<<mu<<' '<<sigma<<' '<<L<<' ';
	for(size_t i=0; i<y.size(); i++)
		out<<y[i]<<' ';
}

/*
#include <fstream>
#include <iomanip>

int main()
{
	RandomNumberGenerator::initialise_instance();

	vector<double> t(100);
	for(size_t i=0; i<t.size(); i++)
		t[i] = pow(1.05, i);

	CAR c(t);
	c.fromPrior();

	fstream fout("output.txt", ios::out);
	fout<<setprecision(12);
	for(int i=0; i<1000; i++)
	{
		CAR c2 = c;
		double logH = c2.perturb();
		if(randomU() <= exp(logH))
			c = c2;
		c.print(fout); fout<<endl;
		cout<<(i+1)<<endl;
	}
	fout.close();

	return 0;
}
*/
