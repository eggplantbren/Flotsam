#include "Limits.h"

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

	isSet = true;
}

