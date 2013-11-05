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
	beta = exp(log(1E-3) + log(1E6)*randomU());
	L = exp(log(L_min) + log(L_max/L_min)*randomU());
	
}

int main()
{
	vector<double> t(3);
	t[0] = 1.; t[1] = 2.; t[2] = 3.;

	CAR c(t);

	return 0;
}

