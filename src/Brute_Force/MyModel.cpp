/*
* Copyright (c) 2009, 2010, 2011, 2012 Brendon J. Brewer.
*
* This file is part of DNest3.
*
* DNest3 is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DNest3 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with DNest3. If not, see <http://www.gnu.org/licenses/>.
*/

#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>

using namespace std;
using namespace DNest3;

MyModel::MyModel()
:tau(Data::get_instance().get_numImages())
,y_qso(Data::get_instance().get_tMin() - 0.5*Data::get_instance().get_tRange(),
	Data::get_instance().get_tMin() + 1.5*Data::get_instance().get_tRange(),
	20000)
,microlensing(Data::get_instance().get_numImages(),
	Curve(Data::get_instance().get_tMin(),
		Data::get_instance().get_tMin() + Data::get_instance().get_tRange(),
		20000))
,noise(Data::get_instance().get_numPoints())
,mu(Data::get_instance().get_numPoints())
{

}

void MyModel::fromPrior()
{
	tau[0] = 0.;
	for(size_t i=1; i<tau.size(); i++)
		tau[i] = 0.1*Data::get_instance().get_tRange()*tan(M_PI*(randomU() - 0.5));

	y_qso.fromPrior();
	for(size_t i=0; i<microlensing.size(); i++)
		microlensing[i].fromPrior();

	noise.fromPrior();

	assemble();
}

double MyModel::perturb()
{
	double logH = 0.;

	int which = randInt(4);

	if(which == 0)
	{
		int which2 = 1 + randInt(tau.size() - 1);
		double u = 0.5 + atan(tau[which2]/(0.1*Data::get_instance().get_tRange()))/M_PI;
		u += pow(10., 1.5 - 6.*randomU())*randn();
		u = mod(u, 1.);
		tau[which2] = 0.1*Data::get_instance().get_tRange()*tan(M_PI*(u - 0.5));
	}
	else if(which == 1)
	{
		int which2 = randInt(microlensing.size());
		logH += microlensing[which2].perturb();
	}
	else if(which == 2)
	{
		logH += y_qso.perturb();
	}
	else if(which == 3)
	{
		logH += noise.perturb();
	}

	assemble();

	return logH;
}

void MyModel::assemble()
{
	const vector<double>& t = Data::get_instance().get_t();
	const vector<int>& id = Data::get_instance().get_ID();

	for(size_t i=0; i<t.size(); i++)
	{
		mu[i] = y_qso.evaluate(t[i] - tau[id[i]]) - y_qso.get_mu()
				+ microlensing[id[i]].evaluate(t[i]);
	}
}



double MyModel::logLikelihood() const
{
	double logL = 0.;

	const vector<double>& y = Data::get_instance().get_y();
	const vector<double>& sig = Data::get_instance().get_sig();

	double s;
	for(size_t i=0; i<y.size(); i++)
	{
		s = noise.get_boost(i)*sig[i];
		logL += -log(s) - 0.5*pow((y[i] - mu[i])/s, 2);
	}

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	for(size_t i=0; i<tau.size(); i++)
		out<<tau[i]<<' ';

	y_qso.print(out); out<<' ';
	for(size_t i=0; i<microlensing.size(); i++)
	{
		microlensing[i].print(out);
		out<<' ';
	}

	noise.print(out);
}

string MyModel::description() const
{
	return string("");
}

