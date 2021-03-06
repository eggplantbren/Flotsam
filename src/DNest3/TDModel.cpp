#include "TDModel.h"
#include "Data.h"
#include "Utils.h"
#include "RandomNumberGenerator.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_linalg.h>

using namespace std;
using namespace DNest3;
using namespace Eigen;

Limits TDModel::limits;

TDModel::TDModel()
:meanVector(Data::get_instance().get_numPoints())
,covarianceMatrix(Data::get_instance().get_numPoints(), Data::get_instance().get_numPoints())
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
	sig_ml.resize(numImages);
	tau_ml.resize(numImages);

	bad_uniforms.resize(numPoints);

	gsl_set_error_handler_off();
}

void TDModel::fromPrior()
{
	for(int i=0; i<numImages; i++)
	{
		if(i==0)
			mag[i] = limits.mag_min + limits.mag_range*randomU();
		else
			mag[i] = tan(M_PI*(randomU() - 0.5));

		if(i == 0)
			tau[i] = 0.;
		else
		{
			// Cauchy prior
			tau[i] = 0.1*limits.tau_range*tan(M_PI*(randomU() - 0.5));
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
	latent_boost1 = randomU();
	boost1 = (latent_boost1 < 0.5)?(1.):(exp(log(1.) + log(100./1.)*(latent_boost1 - 0.5)*2));
	boost2 = exp(log(1.) + log(100./1.)*randomU());

	formMeanVector();
	formCovarianceMatrix();
}

double TDModel::perturb1()
{
	int which = randInt(numImages);

	if(which == 0)
	{
		mag[which] += limits.mag_range
				*pow(10., 1.5 - 6.*randomU())*randn();
		mag[which] = mod(mag[which] - limits.mag_min
					,limits.mag_range)
					+ limits.mag_min;
	}
	else
	{
		mag[which] = atan(mag[which])/M_PI + 0.5;
		mag[which] += pow(10., 1.5 - 6.*randomU())*randn();
		mag[which] = mod(mag[which], 1.);
		mag[which] = tan(M_PI*(mag[which] - 0.5));
	}

	return 0.;
}

double TDModel::perturb2()
{
	double logH = 0.;
	int which = 1 + randInt(numImages - 1);

	tau[which] = atan(tau[which]/(0.1*limits.tau_range))/M_PI + 0.5;
	tau[which] += pow(10., 1.5 - 6.*randomU())*randn();
	tau[which] = mod(tau[which], 1.);
	tau[which] = 0.1*limits.tau_range*tan(M_PI*(tau[which] - 0.5));

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
	int which = randInt(2);
	if(which == 0)
	{
		f_bad += pow(10., 1.5 - 6.*randomU())*randn();
		f_bad = mod(f_bad, 1.);
	}
	else
	{
		latent_boost1 += pow(10., 1.5 - 6.*randomU())*randn();
		latent_boost1 = mod(latent_boost1, 1.);
		boost1 = (latent_boost1 < 0.5)?(1.):(exp(log(1.) + log(100./1.)*(latent_boost1 - 0.5)*2));

		boost2 = log(boost2);
		boost2 += log(100./1.)*pow(10., 1.5 - 6.*randomU())*randn();
		boost2 = mod(boost2 - log(1.), log(100.)) + log(1.);
		boost2 = exp(boost2);
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
	sig_qso = exp(logSig_qso);
	tau_qso = exp(logTau_qso);
	for(int i=0; i<numImages; i++)
	{
		sig_ml[i] = exp(logSig_ml[i]);
		tau_ml[i] = exp(logTau_ml[i]);
	}
	coeff = pow(sig_qso, 2)*tau_qso;

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
		sig *= boost1;
		if(bad_uniforms[i] <= f_bad)
			sig *= boost2;
		covarianceMatrix(i, i) += pow(sig, 2);
	}

	cholesky = covarianceMatrix.llt();
}

void TDModel::formMeanVector()
{
	for(int i=0; i<numPoints; i++)
	{
		if(Data::get_instance().get_ID()[i] == 0)
			meanVector(i) = mag[0];
		else
			meanVector(i) = mag[0] + mag[Data::get_instance().get_ID()[i]];
	}
}

double TDModel::covariance(double t1, double t2, int ID1, int ID2)
{
	double C = coeff*exp(-abs((t1 - tau[ID1]) - (t2 - tau[ID2]))/tau_qso);

	if(ID1 == ID2)
		C += pow(sig_ml[ID1], 2)*exp(-pow(abs(t1 - t2)/tau_ml[ID1], alpha));

	return C;
}

double TDModel::logLikelihood() const
{
	VectorXd y(numPoints);
	for(int i=0; i<numPoints; i++)
		y(i) = Data::get_instance().get_y()[i] - meanVector(i);

	MatrixXd L = cholesky.matrixL();
	double logDeterminant = 0.;
	for(int i=0; i<numPoints; i++)
		logDeterminant += 2.*log(L(i,i));

	// C^-1*(y-mu)
	VectorXd solution = cholesky.solve(y);

	// y*solution
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
	out<<f_bad<<' '<<boost1<<' '<<boost2<<' ';

	for(int i=0; i<numPoints; i++)
		out<<bad_uniforms[i]<<' ';
}

void TDModel::print2(ostream& out) const
{
	out<<setprecision(12);
	// Write out mean vector and covariance matrix
	for(int i=0; i<numPoints; i++)
		out<<meanVector(i)<<' ';
	for(int i=0; i<numPoints; i++)
		for(int j=0; j<numPoints; j++)
			out<<covarianceMatrix(i, j)<<' ';
}

istream& TDModel::read(istream& in)
{
	for(int i=0; i<numImages; i++)
		in>>mag[i];
	for(int i=0; i<numImages; i++)
		in>>tau[i];
	for(int i=0; i<numImages; i++)
		in>>logSig_ml[i];
	for(int i=0; i<numImages; i++)
		in>>logTau_ml[i];

	in>>alpha;
	in>>logSig_qso;
	in>>logTau_qso;
	in>>f_bad;
	in>>boost1;
	in>>boost2;

	for(int i=0; i<numPoints; i++)
		in>>bad_uniforms[i];

	formMeanVector();
	formCovarianceMatrix();

	return in;
}

string TDModel::description() const
{
	return string("mag, tau, logSig_ml, logTau_ml, alpha")
			+ string(", logSig_qso, logTau_qso")
			+ string(", f_bad, boost1, boost2");
}

