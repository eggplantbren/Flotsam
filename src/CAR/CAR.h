#ifndef _CAR_
#define _CAR_

#include <vector>

class CAR
{
	private:
		// Timestamps
		std::vector<double> t;

		// Latent variables with N(0,1) priors
		std::vector<double> n;

		// Hyperparameters
		double mu, beta, L;

	public:
		CAR(const std::vector<double>& times);

		void fromPrior();

};

#endif

