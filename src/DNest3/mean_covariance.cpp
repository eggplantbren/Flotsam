#include <iostream>
#include <fstream>
#include "CommandLineOptions.h"
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

	TDModel t;

	// Load posterior samples
	// and save the mean and covariance to a file
	fstream fin("posterior_sample.txt", ios::in);
	fstream fout("mean_covariance.txt", ios::out);
	if(!fin)
		cerr<<"# ERROR. Couldn't load posterior_sample.txt"<<endl;
	int k = 0;
	while(t.read(fin))
	{
		t.print2(fout); fout<<endl;
		cout<<++k<<endl;
		if(k >= 16)
			break;
	}
	fout.close();
	fin.close();

	return 0;
}

