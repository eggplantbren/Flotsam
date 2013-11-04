#ifndef _Curve_
#define _Curve_

#include <vector>
#include <ostream>

class Curve
{
	private:
		int N;
		double t_min, t_max, t_range, dt;
		std::vector<double> n, y;

		// Mean, short-term sd and timescale
		double mu, beta, L;

		double log_prob() const;
		void assemble();
		void disassemble();
		double perturb1();
		double perturb2();
		double perturb_param();

	public:
		Curve(double t_min, double t_max, int N);

		void fromPrior();
		double perturb();
		void print(std::ostream& out) const;

		double evaluate(double t) const;
};

#endif

