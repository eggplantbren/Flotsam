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
:light_curve(Data::get_instance().get_t())
{

}

void MyModel::fromPrior()
{
	light_curve.fromPrior();
}

double MyModel::perturb()
{
	double logH = light_curve.perturb();
	return logH;
}

double MyModel::logLikelihood() const
{
	const vector<double>& y = light_curve.get_y();
	const vector<double>& Y = Data::get_instance().get_y();
	const vector<double>& sig = Data::get_instance().get_sig();

	double logL = 0.;

	for(size_t i=0; i<y.size(); i++)
		logL += -log(sig[i]) - 0.5*pow((y[i] - Y[i])/sig[i], 2);

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	light_curve.print(out);
}

string MyModel::description() const
{
	return string("");
}

