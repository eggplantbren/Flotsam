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
	alpha += limits.alpha_range*pow(10., 1.5 - 6.*randomU())*randn();
	alpha = mod(alpha - limits.alpha_min, limits.alpha_range)
			+ limits.alpha_min;
	return 0.;
}

double TDModel::perturb6()
{
	logSig_qso += limits.logSig_qso_range
				*pow(10., 1.5 - 6.*randomU())*randn();
	logSig_qso = mod(logSig_qso - limits.logSig_qso_min,
				limits.logSig_qso_range) + limits.logSig_qso_min;
	return 0.;
}

double TDModel::perturb7()
{
	logTau_qso += limits.logTau_qso_range
				*pow(10., 1.5 - 6.*randomU())*randn();
	logTau_qso = mod(logTau_qso - limits.logTau_qso_min,
				limits.logTau_qso_range) + limits.logTau_qso_min;
	return 0.;
}

double TDModel::perturb()
{
	double logH = 0.;

	int which = randInt(7);
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

	// formMeanVector();
	formCovarianceMatrix();
	return logH;
}

void TDModel::formCovarianceMatrix()
{
	covarianceMatrix.setZero();

	// Fill covariance matrix
	// with covariance function evaluations
	for(int i=0; i<Data::get_instance().get_numPoints(); i++)
	{
		for(int j=i; j<Data::get_instance().get_numPoints(); j++)
		{
			covarianceMatrix(i, j) =
				covariance(Data::get_instance().get_t()[i],
					   Data::get_instance().get_t()[j],
					   Data::get_instance().get_ID()[i],
					   Data::get_instance().get_ID()[j]);

			if(i != j)
				covarianceMatrix(j, i) = covarianceMatrix(i, j);
		}
	}

	// Add diagonal noise component
	for(int i=0; i<Data::get_instance().get_numPoints(); i++)
		covarianceMatrix(i, i) += pow(Data::get_instance().get_sig()[i], 2);

	cholesky = covarianceMatrix.llt();
}


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

double TDModel::logLikelihood() const
{
	Matrix L = cholesky.matrixL();

	Vector y(Data::get_instance().get_numPoints());
	for(int i=0; i<Data::get_instance().get_numPoints(); i++)
	{
		y(i) = Data::get_instance().get_y()[i]
			-mag[Data::get_instance().get_ID()[i]];
	}

	double logDeterminant = 0.;
	for(int i=0; i<Data::get_instance().get_numPoints(); i++)
		logDeterminant += 2.*log(L(i,i));

	Vector solution = cholesky.solve(y); // C^-1*(y-mu)
	double exponent = y.dot(solution);

	return -0.5*Data::get_instance().get_numPoints()*log(2*M_PI)
			- 0.5*logDeterminant - 0.5*exponent;
}

