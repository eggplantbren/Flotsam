#ifndef _TDModel_
#define _TDModel_

#include "Model.h"
#include "Data.h"
#include "Limits.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

class TDModel:public DNest3::Model
{
	private:
		// Prior bounds on everything
		static Limits limits;

		// Useful integers
		int numPoints, numImages;

		// Mean magnitude of image 0, followed by magnitudes of
		// other images relative to image 0
		std::vector<double> mag;

		// Time delays (first is zero by definition)
		std::vector<double> tau;		

		// Microlensing amplitudes
		std::vector<double> logSig_ml;

		// Microlensing timescales
		std::vector<double> logTau_ml;

		// Microlensing smoothness
		double alpha;

		// QSO Variability amplitude IN UNITS OF SQRT(TAU)
		double logSig_qso;

		// QSO variability timescale
		double logTau_qso;

		// Hyperparameters for error bar boosts
		double f_bad;
		double latent_boost1, boost1, boost2;
		std::vector<double> bad_uniforms; // U(0,1) prior


		// Temporary stuff
		std::vector<double> sig_ml, tau_ml;
		double sig_qso, tau_qso, coeff;

		// Mean vector
		Eigen::VectorXd meanVector;

		// Covariance matrix and its Cholesky decomposition
		Eigen::MatrixXd covarianceMatrix;
		Eigen::LLT<Eigen::MatrixXd> cholesky;

		// Noise-free covariance function
		double covariance(double t1, double t2, int ID1, int ID2);

		// Helper methods
		double perturb1();
		double perturb2();
		double perturb3();
		double perturb4();
		double perturb5();
		double perturb6();
		double perturb7();
		double perturb8();
		double perturb9();

		void formMeanVector();
		void formCovarianceMatrix();

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

		// Print mean vector and covariance matrix
		void print2(std::ostream& out) const;

		// Read from stream
		std::istream& read(std::istream& in);

		// Return string with column information
		std::string description() const;
};

#endif

