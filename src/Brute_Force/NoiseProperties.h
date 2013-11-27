#ifndef _NoiseProperties_
#define _NoiseProperties_

#include <vector>
#include <ostream>

class NoiseProperties
{
	private:
		double nu;
		double latent_boost;

	public:
		NoiseProperties();

		void fromPrior();
		double perturb();

		double get_nu() const { return nu; }
		double get_boost() const;
		void print(std::ostream& out) const;
};

#endif

