#ifndef _NoiseProperties_
#define _NoiseProperties_

#include <vector>

class NoiseProperties
{
	private:
		std::vector<double> u;
		double f_bad, boost1, boost2;

	public:
		NoiseProperties(int N);

		void fromPrior();
		double perturb();

		double get_boost(int i) const;
};

#endif

