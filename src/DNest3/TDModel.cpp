#include "TDModel.h"
#include "Data.h"
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
	}

	alpha = limits.alpha_min + limits.alpha_range*randomU();
	logSig_qso = limits.logSig_qso_min + limits.logSig_qso_range*randomU();
	logTau_qso = limits.logTau_qso_min + limits.logTau_qso_range*randomU();
}


/*double TDModel::perturbHelper1()
{
	int which = 1+randInt(timeDelays.size()-1);
	timeDelays[which] += data.tRange*pow(10.0, 1.5-6*randomU())*randn();
	timeDelays[which] = mod(timeDelays[which] + 0.5*data.tRange, data.tRange) - 0.5*data.tRange;
	return 0;
}

double TDModel::perturbHelper2()
{
	int which = randInt(sigMicrolensing.size());
	double temp = log(sigMicrolensing[which]);
	temp += log(maxSig/minSig)*pow(10.0, 1.5-6*randomU())*randn();
	temp = mod(temp - log(minSig), log(maxSig/minSig)) + log(minSig);
	sigMicrolensing[which] = exp(temp);
	return 0.0;
}

double TDModel::perturbHelper3()
{
	int which = randInt(tauMicrolensing.size());
	double temp = log(tauMicrolensing[which]);
	temp += log(1E4)*pow(10.0, 1.5-6*randomU())*randn();
	temp = mod(temp - log(1E-2*data.tRange), log(1E4)) + log(1E-2*data.tRange);
	tauMicrolensing[which] = exp(temp);
	return 0.0;
}

double TDModel::perturbHelper4()
{
	int which = randInt(meanMagnitudes.size());
	meanMagnitudes[which] += data.yRange*pow(10.0, 1.5-6*randomU())*randn();
	meanMagnitudes[which] = mod(meanMagnitudes[which] - data.yMin, data.yRange) + data.yMin;
	return 0.0;
}

double TDModel::perturbHelper5()
{
	double temp = log(sigIntrinsic);
	temp += log(maxSig/minSig)*pow(10.0, 1.5-6*randomU())*randn();
	temp = mod(temp - log(minSig), log(maxSig/minSig)) + log(minSig);
	sigIntrinsic = exp(temp);
	return 0.0;	
}

double TDModel::perturbHelper6()
{
	double temp = log(tauIntrinsic);
	temp += log(1E4)*pow(10.0, 1.5-6*randomU())*randn();
	temp = mod(temp - log(1E-2*data.tRange), log(1E4)) + log(1E-2*data.tRange);
	tauIntrinsic = exp(temp);
	return 0.0;
}

double TDModel::perturbHelper7()
{
	alphaMicrolensing += pow(10.0, 1.5-6*randomU())*randn();
	alphaMicrolensing = mod(alphaMicrolensing-1.0, 1.0) + 1.0;
	return 0.0;
}

double TDModel::perturbHelper8()
{
	sigmaBoost = log(sigmaBoost);
	sigmaBoost += log(10.0)*pow(10.0, 1.5-6*randomU())*randn();
	sigmaBoost = mod(sigmaBoost, log(10.0));
	sigmaBoost = exp(sigmaBoost);
	return 0.0;
}

double TDModel::perturbHelper9()
{
	alphaIntrinsic += pow(10.0, 1.5-6*randomU())*randn();
	alphaIntrinsic = mod(alphaIntrinsic-1.0, 1.0) + 1.0;
	return 0.0;
}


double TDModel::perturb()
{
	int which = (randomU() < 0.3)?(0):randInt(9);
	double logh = 0;
	bool necessary = true;

	switch(which)
	{
		case 0:
			logh += perturbHelper1();
			break;
		case 1:
			logh += perturbHelper2();
			break;
		case 2:
			logh += perturbHelper3();
			break;
		case 3:
			logh += perturbHelper4();
			necessary = false;
			break;
		case 4:
			logh += perturbHelper5();
			break;
		case 5:
			logh += perturbHelper6();
			break;
		case 6:
			logh += perturbHelper7();
			break;
		case 7:
			logh += perturbHelper8();
			break;
		case 8:
			logh += perturbHelper9();
			break;
	}

	DNest::Model::perturb();
	if(necessary)
		formCovarianceMatrix();
	calculateLogLikelihood();
	return logh;
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

double TDModel::covariance(double t1, double t2, int qsoID1, int qsoID2)
{
	double result = pow(sigIntrinsic, 2)*exp(-pow(abs((t1 - timeDelays[qsoID1]) - (t2 - timeDelays[qsoID2]))/tauIntrinsic, alphaIntrinsic));
	if(qsoID1 == qsoID2)
		result += pow(sigMicrolensing[qsoID1], 2)*exp(-pow(abs(t1 - t2)/tauMicrolensing[qsoID1], alphaMicrolensing));
	return result;
}

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

