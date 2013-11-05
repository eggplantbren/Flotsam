#include "CAR.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include "RandomNumberGenerator.h"

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

	for(size_t i=0; i<y.size(); i++)
		cout<<y[i]<<endl;
}

int main()
{
	RandomNumberGenerator::initialise_instance();

	vector<double> t(100);
	for(size_t i=0; i<t.size(); i++)
		t[i] = i;
	t[98] = 105;
	t[99] = 110;

	CAR c(t);
	c.fromPrior();

	return 0;
}

