#include <iostream>
#include "Start.h"
#include "Data.h"
#include "TDModel.h"

using namespace std;
using namespace DNest3;

int main(int argc, char** argv)
{
	// Process command line options and load data
	CommandLineOptions options(argc, argv);
	string dataFile = options.get_dataFile();
	Data::get_instance().load(dataFile.c_str());

	// Initialise the sampler
	MTSampler<TDModel> sampler = setup_mt<TDModel>(options);

	// Go!
	sampler.run();

	return 0;
}

