#ifndef _Limits_
#define _Limits_

#include "Data.h"

/*
* Define prior bounds on the parameters.
*/

class Limits
{
	friend class TDModel;

	private:
		// Naming convention: varName_min, varName_max, varName_range

		// Mean magnitudes
		std::vector<double> mag_min, mag_max, mag_range;

		// Time delays (first is zero by definition)
		std::vector<double> tau_min, tau_max, tau_range;

		// Microlensing amplitudes
		std::vector<double> logSig_ml_min, logSig_ml_max,
						logSig_ml_range;

		// Microlensing smoothness
		double alpha_min, alpha_max, alpha_range;

		// QSO Variability amplitude
		double logSig_qso_min, logSig_qso_max, logSig_qso_range;

		// QSO variability timescale
		double logTau_qso_min, logTau_qso_max, logTau_qso_range;
		

	public:
		// Default constructor: Do nothing
		Limits();

		// Use data to construct limits
		void set(const Data& data);

};

#endif

