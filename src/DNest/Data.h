#ifndef _Data_
#define _Data_

#include <vector>

class Data
{
	private:
		// The data
		std::vector<double> t, y, sig;
		std::vector<int> ID;

		bool loaded;

		// Some summary statistics
		int numPoints, numImages;
		double tMin, tMax, tRange;
		std::vector<double> yMean, yStDev; // One for each image

		void computeSummaries();

	public:
		Data();
		void load(const char* filename);

		// A bunch of getters
		double get_tMin() const { return tMin; }
		double get_tRange() const { return tRange; }
		int get_numPoints() const { return numPoints; }
		int get_numImages() const { return numImages; }

};

#endif

