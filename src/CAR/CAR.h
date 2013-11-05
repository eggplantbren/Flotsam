#ifndef _CAR_
#define _CAR_

#include <vector>

class CAR
{
	private:
		// Latent variables with N(0,1) priors
		std::vector<double> n;

	public:
		CAR(int N);

};

#endif

