#include "TDModel.h"
#include "Data.h"
#include "Utils.h"
#include "RandomNumberGenerator.h"
#include <iostream>
#include <cmath>
#include <gsl/gsl_linalg.h>

using namespace std;
using namespace DNest3;

Limits TDModel::limits;

TDModel::TDModel()
:meanVector(1)
,covarianceMatrix(1, 1)
,cholesky(1, 1)
{
	if(!Data::get_instance().get_loaded())
	{
		cerr<<"# Data has not been loaded! Cannot construct TDModel."<<endl;
		exit(0);
	}
	limits.set(Data::get_instance());

	numPoints = Data::get_instance().get_numPoints();
	numImages = Data::get_instance().get_numImages();

	mag.resize(numImages);
	tau.resize(numImages);
	logSig_ml.resize(numImages);
	logTau_ml.resize(numImages);

	meanVector = Vector(numPoints);
	covarianceMatrix = Matrix(numPoints, numPoints);
	cholesky = Matrix(numPoints, numPoints);
	bad_uniforms.resize(numPoints);

	gsl_set_error_handler_off();
}

void TDModel::fromPrior()
{
	for(int i=0; i<numImages; i++)
	{
		mag[i] = limits.mag_min[i] + limits.mag_range[i]*randomU();

		if(i == 0)
			tau[i] = 0.;
		else
		{
			// Cauchy prior with bounds
			do
			{
				tau[i] = limits.tau_min +
						limits.tau_range*randomU();
			}while(randomU() >= 1./(1. + pow(tau[i]/
						(0.1*limits.tau_range), 2)));
		}

		logSig_ml[i] = limits.logSig_ml_min[i]
					+ limits.logSig_ml_range[i]*randomU();
		logTau_ml[i] = limits.logTau_ml_min[i]
					+ limits.logTau_ml_range[i]*randomU();
	}

	alpha = limits.alpha_min + limits.alpha_range*randomU();
	logSig_qso = limits.logSig_qso_min + limits.logSig_qso_range*randomU();
	logTau_qso = limits.logTau_qso_min + limits.logTau_qso_range*randomU();


	for(int i=0; i<numPoints; i++)
		bad_uniforms[i] = randomU();
	f_bad = randomU();
	boost = exp(log(1.) + log(100./1.)*randomU());

	formMeanVector();
	formCovarianceMatrix();
}

double TDModel::perturb1()
{
	int which = randInt(numImages);

	mag[which] += limits.mag_range[which]
			*pow(10., 1.5 - 6.*randomU())*randn();
	mag[which] = mod(mag[which] - limits.mag_min[which]
				,limits.mag_range[which])
				+ limits.mag_min[which];
	return 0.;
}

double TDModel::perturb2()
{
	double logH = 0.;
	int which = 1 + randInt(numImages - 1);

	logH -= -log(1. + pow(tau[which]/(0.1*limits.tau_range), 2));
	tau[which] += limits.tau_range
			*pow(10., 1.5 - 6.*randomU())*randn();
	tau[which] = mod(tau[which] - limits.tau_min
				,limits.tau_range)
				+ limits.tau_min;
	logH += -log(1. + pow(tau[which]/(0.1*limits.tau_range), 2));

	return logH;
}

double TDModel::perturb3()
{
	int which = randInt(numImages);
	logSig_ml[which] += limits.logSig_ml_range[which]
			*pow(10., 1.5 - 6.*randomU())*randn();
	logSig_ml[which] = mod(logSig_ml[which] - limits.logSig_ml_min[which]
				,limits.logSig_ml_range[which])
				+ limits.logSig_ml_min[which];
	return 0.;
}

double TDModel::perturb4()
{
	int which = randInt(numImages);
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

double TDModel::perturb8()
{
	double scale = pow(10., 1.5 - 6.*randomU());
	int which = randInt(2);
	if(which == 0)
	{
		f_bad += scale*randn();
		f_bad = mod(f_bad, 1.);
	}
	else
	{
		boost = log(boost);
		boost += log(100./1.)*scale*randn();
		boost = mod(boost - log(1.), log(100.)) + log(1.);
		boost = exp(boost);
	}

	return 0.;
}

double TDModel::perturb9()
{
	double logH = 0.;

	// Resample some uniforms
	double chance = pow(10., 0.5 - 4.*randomU());
	for(int i=0; i<numPoints; i++)
	{
		if(randomU() <= chance)
			bad_uniforms[i] = randomU();
	}

	return logH;
}

double TDModel::perturb()
{
	double logH = 0.;

	// Flag -- whether to do each kind of proposal
	vector<bool> change(9, false);
	int num = 0;
	double chance = pow(10., 0.5 - 3.*randomU());
	for(size_t i=0; i<change.size(); i++)
	{
		if(randomU() <= chance)
		{
			change[i] = true;
			num++;
		}
	}
	if(num == 0)
	{
		change[randInt(change.size())] = true;
		num++; 
	}

	if(change[0])
		logH += perturb1();
	if(change[1])
		logH += perturb2();
	if(change[2])
		logH += perturb3();
	if(change[3])
		logH += perturb4();
	if(change[4])
		logH += perturb5();
	if(change[5])
		logH += perturb6();
	if(change[6])
		logH += perturb7();
	if(change[7])
		logH += perturb8();
	if(change[8])
		logH += perturb9();


	formMeanVector();
	formCovarianceMatrix();
	return logH;
}

void TDModel::formCovarianceMatrix()
{
	// Fill covariance matrix
	// with covariance function evaluations
	for(int i=0; i<numPoints; i++)
	{
		for(int j=i; j<numPoints; j++)
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
	double sig;
	for(int i=0; i<numPoints; i++)
	{
		sig = Data::get_instance().get_sig()[i];
		if(bad_uniforms[i] <= f_bad)
			sig *= boost;
		covarianceMatrix(i, i) += pow(sig, 2);
	}

	cholesky = covarianceMatrix;
	gsl_linalg_cholesky_decomp(cholesky.get_gsl_matrix());
}

void TDModel::formMeanVector()
{
	for(int i=0; i<numPoints; i++)
		meanVector(i) = mag[Data::get_instance().get_ID()[i]];
}

double TDModel::covariance(double t1, double t2, int ID1, int ID2)
{
	double sig_qso = exp(logSig_qso);
	double tau_qso = exp(logTau_qso);

	double exponent = abs((t1 - tau[ID1]) - (t2 - tau[ID2]))/tau_qso;
	double C = pow(sig_qso, 2)
				*exp(-exponent);

	vector<double> sig_ml(numImages);
	vector<double> tau_ml(numImages);
	for(int i=0; i<numImages; i++)
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
	Vector y(numPoints);
	for(int i=0; i<numPoints; i++)
		y(i) = Data::get_instance().get_y()[i] - meanVector(i);

	double logDeterminant = 0.;
	for(int i=0; i<numPoints; i++)
		logDeterminant += 2.*log(cholesky(i,i));

	// C^-1*(y-mu)
	Vector solution(numPoints);
	gsl_linalg_cholesky_solve(cholesky.get_gsl_matrix(), y.get_gsl_vector(), solution.get_gsl_vector());

	// y . solution
	double exponent = 0.;
	for(int i=0; i<numPoints; i++)
		exponent += y(i)*solution(i);

	double logL = -0.5*numPoints*log(2*M_PI)
			- 0.5*logDeterminant - 0.5*exponent;

	if(isnan(logL) || isinf(logL))
		logL = -1E300;

	return logL;
}

void TDModel::print(ostream& out) const
{
	for(int i=0; i<numImages; i++)
		out<<mag[i]<<' ';
	for(int i=0; i<numImages; i++)
		out<<tau[i]<<' ';
	for(int i=0; i<numImages; i++)
		out<<logSig_ml[i]<<' ';
	for(int i=0; i<numImages; i++)
		out<<logTau_ml[i]<<' ';

	out<<alpha<<' ';
	out<<logSig_qso<<' ';
	out<<logTau_qso<<' ';
	out<<f_bad<<' '<<boost<<' ';
}

string TDModel::description() const
{
	return string("mag, tau, logSig_ml, logTau_ml, alpha")
			+ string(", logSig_qso, logTau_qso")
			+ string(", f_bad, boost");
}

