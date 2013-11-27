#ifndef _CAR_
#define _CAR_

#include "ArgSorter.h"
#include <vector>
#include <ostream>

class CAR
{
	private:
		// Allow mean to be free or fixed?
		bool zero_mean;

		// Timestamps
		ArgSorter<double> t;

		// Latent variables with N(0,1) priors
		std::vector<double> n;

		// Actual value of function
		std::vector<double> y;

		// Hyperparameters
		double mu, beta, L;

		// Limits on L
		double L_min, L_max;

		// Compute y from n (and hyperparameters)
		void assemble();

	public:
		CAR(bool zero_mean, const std::vector<double>& times);

		void fromPrior();
		double perturb();

		void print(std::ostream& out) const;

		const std::vector<double>& get_y() const { return y; }
};

#endif

