#include "TDModel.h"
#include "Data.h"
#include "Utils.h"
#include "RandomNumberGenerator.h"
#include <iostream>

using namespace std;
using namespace DNest3;

TDModel::TDModel()
{
	if(!Data::get_instance().get_loaded())
	{
		cerr<<"# Data has not been loaded! Cannot construct TDModel"<<endl;
		exit(0);
	}
}

void TDModel::fromPrior()
{
	for(int i=0; i<Data::get_instance().get_numImages(); i++)
	{
		mag[i] = limits.mag_min[i] + limits.mag_range[i]*randomU();
		tau[i] = limits.tau_min[i] + limits.tau_range[i]*randomU();
		logSig_ml[i] = limits.logSig_ml_min[i]
					+ limits.logSig_ml_range[i]*randomU();
		logTau_ml[i] = limits.logTau_ml_min[i]
					+ limits.logTau_ml_range[i]*randomU();
	}

	alpha = limits.alpha_min + limits.alpha_range*randomU();
	logSig_qso = limits.logSig_qso_min + limits.logSig_qso_range*randomU();
	logTau_qso = limits.logTau_qso_min + limits.logTau_qso_range*randomU();
}

double TDModel::perturb1()
{
	int which = randInt(Data::get_instance().get_numImages());
	mag[which] += limits.mag_range[which]
			*pow(10., 1.5 - 6.*randomU())*randn();
	mag[which] = mod(mag[which] - limits.mag_min[which]
				,limits.mag_range[which])
				+ limits.mag_min[which];
	return 0.;
}

double TDModel::perturb2()
{
	int which = randInt(Data::get_instance().get_numImages());
	tau[which] += limits.tau_range[which]
			*pow(10., 1.5 - 6.*randomU())*randn();
	tau[which] = mod(tau[which] - limits.tau_min[which]
				,limits.tau_range[which])
				+ limits.tau_min[which];
	return 0.;
}

double TDModel::perturb3()
{
	int which = randInt(Data::get_instance().get_numImages());
	logSig_ml[which] += limits.logSig_ml_range[which]
			*pow(10., 1.5 - 6.*randomU())*randn();
	logSig_ml[which] = mod(logSig_ml[which] - limits.logSig_ml_min[which]
				,limits.logSig_ml_range[which])
				+ limits.logSig_ml_min[which];
	return 0.;
}

double TDModel::perturb4()
{
	int which = randInt(Data::get_instance().get_numImages());
	logTau_ml[which] += limits.logTau_ml_range[which]
			*pow(10., 1.5 - 6.*randomU())*randn();
	logTau_ml[which] = mod(logTau_ml[which] - limits.logTau_ml_min[which]
				,limits.logTau_ml_range[which])
				+ limits.logTau_ml_min[which];
	return 0.;
}

double TDModel::perturb5()
{
	double alpha += limits.alpha_range*pow(10., 1.5 - 6.*randomU())*randn();
	alpha = mod(alpha - limits.alpha_min, limits.alpha_range)
			+ limits.alpha_min;
	return 0.;
}

double TDModel::perturb6()
{
	double logSig_qso += limits.logSig_qso_range
				*pow(10., 1.5 - 6.*randomU())*randn();
	logSig_qso = mod(logSig_qso - limits.logSig_qso_min,
				limits.logSig_qso_range) + limits.logSig_qso_min;
	return 0.;
}

double TDModel::perturb7()
{
	double logTau_qso += limits.logTau_qso_range
				*pow(10., 1.5 - 6.*randomU())*randn();
	logTau_qso = mod(logTau_qso - limits.logTau_qso_min,
				limits.logTau_qso_range) + limits.logTau_qso_min;
	return 0.;
}

double TDModel::perturb()
{
	int which = randInt(7);
	double logH = 0;
	bool necessary = true;

	switch(which)
	{
		case 0:
			logH += perturb1();
			break;
		case 1:
			logH += perturb2();
			break;
		case 2:
			logH += perturb3();
			break;
		case 3:
			logH += perturb4();
			break;
		case 4:
			logH += perturb5();
			break;
		case 5:
			logH += perturb6();
			break;
		case 6:
			logH += perturb7();
			break;
	}

	return logH;
}

void TDModel::formCovarianceMatrix()
{
	covMat.setZero();

	// Fill covariance matrix
	// with covariance function evaluations
	for(int i=0; i<data.N; i++)
	{
		for(int j=i; j<data.N; j++)
		{
			covMat(i, j) = covariance(data.t[i], data.t[j], data.qsoID[i], data.qsoID[j]);
			if(i != j)
				covMat(j, i) = covMat(i, j);
		}
	}

	// Add diagonal noise component
	for(int i=0; i<data.N; i++)
		covMat(i,i) += pow(data.sig[i]*sigmaBoost, 2);

	cholesky = covMat.llt();
}
*/
double TDModel::covariance(double t1, double t2, int ID1, int ID2)
{
	double sig_qso = exp(logSig_qso);
	double tau_qso = exp(logTau_qso);

	double exponent = abs((t1 - tau[ID1]) - (t2 - tau[ID2]))/tau_qso;
	double C = pow(sig_qso, 2)
				*exp(-exponent);

	vector<double> sig_ml(Data::get_instance().get_numImages());
	vector<double> tau_ml(Data::get_instance().get_numImages());
	for(int i=0; i<Data::get_instance().get_numImages(); i++)
	{
		sig_ml[i] = exp(logSig_ml[i]);
		tau_ml[i] = exp(logTau_ml[i]);
	}

	if(ID1 == ID2)
	{
		exponent = pow(abs(t1 - t2)/tau_ml[ID1], alpha);
		C += pow(sig_ml[ID1], 2)*exp(-exponent);
	}
	return C;
}
/*
void TDModel::calculateLogLikelihood()
{
	MatrixXd L = cholesky.matrixL();

	VectorXd y(data.N);
	for(int i=0; i<data.N; i++)
		y(i) = data.y[i] - meanMagnitudes[data.qsoID[i]];

	double logDeterminant = 0;
	for(int i=0; i<data.N; i++)
		logDeterminant += 2*log(L(i,i));

	VectorXd solution = cholesky.solve(y); // C^-1*(y-mu)
	double exponent = y.dot(solution);

	double loglGood = -0.5*data.N*log(2*pi) - 0.5*logDeterminant - 0.5*exponent;

	logl.logl = loglGood;
}

double TDModel::getLogLikelihood() const
{
	return logl.logl;
}

void TDModel::print(ostream& out) const
{
	for(size_t i=0; i<timeDelays.size(); i++)
		out<<timeDelays[i]<<' ';
	for(size_t i=0; i<sigMicrolensing.size(); i++)
		out<<sigMicrolensing[i]<<' ';
	for(size_t i=0; i<tauMicrolensing.size(); i++)
		out<<tauMicrolensing[i]<<' ';
	for(size_t i=0; i<meanMagnitudes.size(); i++)
		out<<meanMagnitudes[i]<<' ';
	out<<sigIntrinsic<<' '<<tauIntrinsic<<' '<<alphaIntrinsic<<' '<<alphaMicrolensing<<' '<<sigmaBoost<<' '<<logl.logl;
}

void TDModel::read(istream& in)
{
	for(size_t i=0; i<timeDelays.size(); i++)
		in>>timeDelays[i];
	for(size_t i=0; i<sigMicrolensing.size(); i++)
		in>>sigMicrolensing[i];
	for(size_t i=0; i<tauMicrolensing.size(); i++)
		in>>tauMicrolensing[i];
	for(size_t i=0; i<meanMagnitudes.size(); i++)
		in>>meanMagnitudes[i];
	in>>sigIntrinsic;
	in>>tauIntrinsic;
	in>>alphaIntrinsic;
	in>>alphaMicrolensing;
	in>>sigmaBoost;
	in>>logl.logl;
}

void TDModel::loadData(const char* filename)
{
	data.load(filename);
}

vector<double> TDModel::evaluate(const vector<double>& tEvaluations, const vector<int>& qsoID)
{
	assert(tEvaluations.size() == qsoID.size());

	// Make covariance matrix of joint prior for curve and data
	unsigned int N = data.t.size();
	vector<double> allTimes(N + tEvaluations.size());
	vector<int> allQsoIDs(N + tEvaluations.size());
	for(size_t i=0; i<N; i++)
	{
		allTimes[i] = data.t[i];
		allQsoIDs[i] = data.qsoID[i];
	}
	for(size_t i=0; i<tEvaluations.size(); i++)
	{
		allTimes[i + N] = tEvaluations[i];
		allQsoIDs[i+N] = qsoID[i];
	}
	
	MatrixXd bigC;
	bigC.setZero(allTimes.size(), allTimes.size());
	for(size_t i=0; i<allTimes.size(); i++)
	{
		for(size_t j=i; j<allTimes.size(); j++)
		{
			bigC(i, j) = covariance(allTimes[i], allTimes[j], allQsoIDs[i], allQsoIDs[j]);
			bigC(j, i) = bigC(i, j);
		}
	}
	for(size_t i=0; i<N; i++)
		bigC(i, i) += pow(data.sig[i]*sigmaBoost, 2);

	// Submatrices
	MatrixXd KXXstar = bigC.block(0, N, N, tEvaluations.size());
	MatrixXd KXstarX = bigC.block(N, 0, tEvaluations.size(), N);
	MatrixXd KXstarXstar = bigC.block(N, N, tEvaluations.size(), tEvaluations.size());

	MatrixXd A = bigC.block(0, 0, N, N).inverse();
	MatrixXd C = KXstarXstar - KXstarX*A*KXXstar;

	vector<double> meanVec1 = data.y;
	for(size_t i=0; i<data.y.size(); i++)
		meanVec1[i] -= meanMagnitudes[data.qsoID[i]];
	vector<double> meanVec2(tEvaluations.size());
	for(size_t i=0; i<meanVec2.size(); i++)
		meanVec2[i] = meanMagnitudes[qsoID[i]];
	VectorXd _meanVec1 = VectorXd::Map(&meanVec1[0], meanVec1.size());
	VectorXd _meanVec2 = VectorXd::Map(&meanVec2[0], meanVec2.size());
	VectorXd meanVec = _meanVec2 + KXstarX*A*_meanVec1;

	vector<double> meanVec_(tEvaluations.size());
	for(size_t i=0; i<tEvaluations.size(); i++)
		meanVec_[i] = meanVec(i);

	// Do Cholesky decomposition
	MatrixXd L = C.llt().matrixL();

	// Make references to things so we can do linear algebra on them
	VectorXd normals(tEvaluations.size());
	for(int i=0; i<normals.size(); i++)
		normals(i) = randn();
	VectorXd y = L*normals + meanVec;

	vector<double> result(tEvaluations.size());
	for(size_t i=0; i<result.size(); i++)
		result[i] = y(i);

	return result;
}

Data& TDModel::getData()
{
	return data;
}
*/

