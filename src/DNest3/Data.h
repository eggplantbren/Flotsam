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

		// Singleton
		static Data instance;

	public:
		Data();
		void load(const char* filename);

		// A bunch of getters
		double get_tMin() const { return tMin; }
		double get_tRange() const { return tRange; }
		int get_numPoints() const { return numPoints; }
		int get_numImages() const { return numImages; }
		const std::vector<double>& get_yMean() const { return yMean; }
		const std::vector<double>& get_yStDev() const { return yStDev; }

		static Data get_instance() { return instance; }

};

#endif

