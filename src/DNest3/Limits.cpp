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

	tau_min.resize(data.get_numImages());
	tau_max.resize(data.get_numImages());
	tau_range.resize(data.get_numImages());
	for(int i=0; i<data.get_numImages(); i++)
	{
		tau_min[i] = -data.get_tRange();
		tau_max[i] =  data.get_tRange();
		tau_range[i] = tau_max[i] - tau_min[i];
	}

	logSig_ml_min.resize(data.get_numImages());
	logSig_ml_max.resize(data.get_numImages());
	logSig_ml_range.resize(data.get_numImages());
	for(int i=0; i<data.get_numImages(); i++)
	{
		logSig_ml_min[i] = log(1E-2*data.get_yStDev()[i]);
		logSig_ml_max[i] = log(1E+2*data.get_yStDev()[i]);
		logSig_ml_range[i] = logSig_ml_max[i] - logSig_ml_min[i];
	}

	isSet = true;
}

