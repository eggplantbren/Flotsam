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
#include <cmath>

using namespace std;
using namespace DNest3;

MyModel::MyModel()
:mu(Data::get_instance().get_numImages())
,qso_light_curve(true, Data::get_instance().get_t())
{

}

void MyModel::fromPrior()
{
	for(size_t i=0; i<mu.size(); i++)
		mu[i] = 10.*tan(M_PI*(randomU() - 0.5));
	qso_light_curve.fromPrior();
}

double MyModel::perturb()
{
	double logH = 0.;

	int which = randInt(2);
	if(which == 0)
		logH += qso_light_curve.perturb();
	else
	{
		int i = randInt(mu.size());

		mu[i] = 0.5 + atan(mu[i]/10.)/M_PI;
		mu[i] += pow(10., 1.5 - 6.*randomU())*randn();
		mu[i] = mod(mu[i], 1.);
		mu[i] = 10.*tan(M_PI*(mu[i] - 0.5));
	}

	return logH;
}

double MyModel::logLikelihood() const
{
	const vector<double>& y = qso_light_curve.get_y();
	const vector<double>& Y = Data::get_instance().get_y();
	const vector<double>& sig = Data::get_instance().get_sig();
	const vector<int>& id = Data::get_instance().get_ID();

	double logL = 0.;

	double m;
	for(size_t i=0; i<y.size(); i++)
	{
		m = mu[0];
		if(id[i] != 0)
			m += mu[id[i]];
		logL += -log(sig[i]) - 0.5*pow((Y[i] - (y[i] + m))/sig[i], 2);
	}

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	for(size_t i=0; i<mu.size(); i++)
		out<<mu[i]<<' ';
	qso_light_curve.print(out);
}

string MyModel::description() const
{
	return string("mu, qso_light_curve");
}

