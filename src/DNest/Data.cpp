#include "Data.h"
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cmath>

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
		cerr<<"# Error opening "<<filename<<"."<<endl;
		exit(0);
	}

	// Skip past comments
	while(fin.peek() == '#')
		fin.ignore(1048576, '\n');

	while(fin>>temp1 && fin>>temp2 && fin>>temp3 && fin>>temp4)
	{
		t.push_back(temp1);
		y.push_back(temp2);
		sig.push_back(temp3);
		ID.push_back(static_cast<int>(temp4));
	}
	cout<<"# Loaded "<<t.size()<<" points from file "<<filename<<"."<<endl;
	fin.close();

	computeSummaries();
	loaded = true;
}

void Data::computeSummaries()
{
	cout<<"# Summary statistics :"<<endl;

	numPoints = t.size();
	numImages = 1 + *max_element(ID.begin(), ID.end());
	tMin = *min_element(t.begin(), t.end());
	tMax = *max_element(t.begin(), t.end());
	tRange = tMax - tMin;
	cout<<"# Number of images = "<<numImages<<endl;
	cout<<"# Time range = "<<tRange<<endl;

	yMean.assign(numImages, 0.);
	vector<double> ySqMean(numImages, 0.);
	yStDev.assign(numImages, 0.);

	for(int k=0; k<numImages; k++)
	{
		double wTot = 0.;
		for(int i=0; i<numPoints; i++)
			if(ID[i] == k)
				wTot += pow(sig[i], -2);

		for(int i=0; i<numPoints; i++)
		{
			if(ID[i] == k)
			{
				yMean[k] += y[i]*pow(sig[i], -2);
				ySqMean[k] += y[i]*y[i]*pow(sig[i], -2);
			}
		}
		yMean[k] /= wTot;
		ySqMean[k] /= wTot;
		yStDev[k] = sqrt(ySqMean[k] - pow(yMean[k], 2));
	}

	cout<<"# Empirical mean levels of images: ";
	for(int k=0; k<numImages; k++)
		cout<<yMean[k]<<' ';
	cout<<"# Empirical standard deviations of images: ";
	for(int k=0; k<numImages; k++)
		cout<<yStDev[k]<<' ';
	
}

