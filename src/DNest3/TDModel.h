#ifndef _TDModel_
#define _TDModel_

#include "Model.h"
#include "Data.h"
#include "Limits.h"
#include <ostream>
#include <vector>
#include <Eigen/Dense>

// Shorter names for Eigen3 types
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::LLT< Eigen::MatrixXd > Cholesky;

class TDModel:public DNest3::Model
{
	private:
		// Prior bounds on everything
		static Limits limits;

		// Mean magnitudes
		std::vector<double> mag;

		// Time delays (first is zero by definition)
		std::vector<double> tau;		

		// Microlensing amplitudes
		std::vector<double> logSig_ml;

		// Microlensing timescales
		std::vector<double> logTau_ml;

		// Microlensing smoothness
		double alpha;

		// QSO Variability amplitude
		double logSig_qso;

		// QSO variability timescale
		double logTau_qso;

		// Mean vector
		Vector meanVector;

		// Covariance matrix and its Cholesky decomposition
		Matrix covarianceMatrix;
		Cholesky cholesky;

		// Noise-free covariance function
		double covariance(double t1, double t2, int ID1, int ID2);

		// Helper methods
		double perturb1();
		double perturb2();

	public:
		TDModel();

		// Generate the point from the prior
		void fromPrior();

		// Metropolis-Hastings proposals
		double perturb();

		// Likelihood function
		double logLikelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif

