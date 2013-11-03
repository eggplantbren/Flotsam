#ifndef _NoiseProperties_
#define _NoiseProperties_

#include <vector>
#include <ostream>

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
		void print(std::ostream& out) const;
};

#endif

