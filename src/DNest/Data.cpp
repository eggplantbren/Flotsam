#include "Data.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>

using namespace std;

Data::Data()
:loaded(false)
{

}


void Data::load(const char* filename)
{
	if(loaded)
		cerr<<"# Already loaded file."<<endl;

	fstream fin(filename, ios::in);
	double temp1, temp2, temp3, temp4;

	if(!fin)
	{	
		cerr<<"Error opening "<<filename<<"."<<endl;
		exit(0);
	}

	// Skip past comments
	while(fin.peek() == '#')
		fin.ignore(1024*1024, '\n');

	while(fin>>temp1 && fin>>temp2 && fin>>temp3 && fin>>temp4)
	{
		t.push_back(temp1);
		y.push_back(temp2);
		sig.push_back(temp3);
		qsoID.push_back((int)temp4);
	}
	cout<<"# Loaded "<<t.size()<<" points from file "<<filename<<"."<<endl;

	fin.close();

	N = t.size();
	numQSOs = 1 + *max_element(qsoID.begin(), qsoID.end());
	tMin = *min_element(t.begin(), t.end());
	tRange = *max_element(t.begin(), t.end()) - tMin;
	yMin = *min_element(y.begin(), y.end());
	yRange = *max_element(y.begin(), y.end()) - yMin;

	loaded = true;
}

double Data::getTMin() const
{
	return tMin;
}

double Data::getTRange() const
{
	return tRange;
}

int Data::getN() const
{
	return N;
}


