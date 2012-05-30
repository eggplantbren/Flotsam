#ifndef _Data_
#define _Data_

#include <vector>

class Data
{
	friend class TDModel;

	private:
		std::vector<double> t, y, sig;
		std::vector<int> qsoID;

		bool loaded;

		int N;
		int numQSOs;
		double tMin, tRange;
		double yMin, yRange;

	public:
		Data();
		void load(const char* filename);
		double getTMin() const;
		double getTRange() const;
		int getN() const;

};

#endif

