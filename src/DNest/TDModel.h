#ifndef _TDModel_
#define _TDModel_

#include "Model.h"
#include "Data.h"
#include <ostream>
#include <istream>
#include <vector>
#include <Eigen/Dense>

class TDModel:public DNest::Model
{
	private:
		std::vector<double> timeDelays;
		std::vector<double> sigMicrolensing;
		std::vector<double> tauMicrolensing;
		std::vector<double> meanMagnitudes;

		double sigIntrinsic, tauIntrinsic, alphaIntrinsic, alphaMicrolensing, sigmaBoost;
		Eigen::MatrixXd covMat;
		Eigen::LLT< Eigen::MatrixXd > cholesky;

		static const double minSig;
		static const double maxSig;

		static Data data;

		void formCovarianceMatrix();

		double perturbHelper1();
		double perturbHelper2();
		double perturbHelper3();
		double perturbHelper4();
		double perturbHelper5();
		double perturbHelper6();
		double perturbHelper7();
		double perturbHelper8();
		double perturbHelper9();

		// Noise-free covariance function
		double covariance(double t1, double t2, int qsoID1, int qsoID2);

	public:
		TDModel();
		~TDModel();
		DNest::Model* factory() const;
		DNest::Model* clone() const;
		void copyFrom(const DNest::Model* other);
		void fromPrior();
		void calculateLogLikelihood();
		double perturb();
		double getValue();
		double getLogLikelihood() const;
		void print(std::ostream& out) const;
		void read(std::istream& in);
		std::vector<double> evaluate(const std::vector<double>& t, const std::vector<int>& qsoID);
		static void loadData(const char* filename);
		static Data& getData();
};

#endif

