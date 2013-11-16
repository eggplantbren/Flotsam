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

		// Actual value of function
		std::vector<double> y;

		// Hyperparameters
		double mu, sigma, L;

		// Limits on L
		double L_min, L_max;

		// Compute y from n (and hyperparameters)
		void assemble();

	public:
		CAR(const std::vector<double>& times);

		void fromPrior();
		double perturb();
};

#endif

