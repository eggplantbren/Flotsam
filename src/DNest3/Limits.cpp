#include "Limits.h"
#include <cmath>

Limits::Limits()
:isSet(false)
{

}

void Limits::set(const Data& data)
{
	mag_min.resize(data.get_numImages());
	mag_max.resize(data.get_numImages());
	mag_range.resize(data.get_numImages());
	for(int i=0; i<data.get_numImages(); i++)
	{
		mag_min[i] = data.get_yMean()[i] - 10.*data.get_yStDev()[i]; 
		mag_max[i] = data.get_yMean()[i] + 10.*data.get_yStDev()[i];
		mag_range[i] = mag_max[i] - mag_min[i];
	}

	tau_min = -data.get_tRange();
	tau_max =  data.get_tRange();
	tau_range = tau_max - tau_min;

	logSig_ml_min.resize(data.get_numImages());
	logSig_ml_max.resize(data.get_numImages());
	logSig_ml_range.resize(data.get_numImages());
	for(int i=0; i<data.get_numImages(); i++)
	{
		logSig_ml_min[i] = log(1E-2*data.get_yStDev()[i]);
		logSig_ml_max[i] = log(1E+2*data.get_yStDev()[i]);
		logSig_ml_range[i] = logSig_ml_max[i] - logSig_ml_min[i];
	}

	logTau_ml_min.resize(data.get_numImages());
	logTau_ml_max.resize(data.get_numImages());
	logTau_ml_range.resize(data.get_numImages());
	for(int i=0; i<data.get_numImages(); i++)
	{
		logTau_ml_min[i] = log(1E-2*data.get_tRange());
		logTau_ml_max[i] = log(1E+3*data.get_tRange());
		logTau_ml_range[i] = logSig_ml_max[i] - logSig_ml_min[i];
	}

	alpha_min = 1.;
	alpha_max = 2.;
	alpha_range = alpha_max - alpha_min;

	logSig_qso_min = 0.;
	logSig_qso_max = 0.;
	for(int i=0; i<data.get_numImages(); i++)
	{
		logSig_qso_min += 1E-2*data.get_yStDev()[i]/data.get_numImages();
		logSig_qso_max += 1E+2*data.get_yStDev()[i]/data.get_numImages();
	}
	logSig_qso_range = logSig_qso_max - logSig_qso_min;

	logTau_qso_min = log(1E-2*data.get_tRange());
	logTau_qso_max = log(1E+3*data.get_tRange());
	logTau_qso_range = logTau_qso_max - logTau_qso_min;

	meanLogSigmaBoost_min = log(1E-1);
	meanLogSigmaBoost_max = log(1E+1);
	meanLogSigmaBoost_range = meanLogSigmaBoost_max
					- meanLogSigmaBoost_min;

	logStdevLogSigmaBoost_min = log(1E-2);
	logStdevLogSigmaBoost_max = log(1E+1);
	logStdevLogSigmaBoost_range = logStdevLogSigmaBoost_max
					- logStdevLogSigmaBoost_min;

	isSet = true;
}

