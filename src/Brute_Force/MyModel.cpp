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
:delta_mag(Data::get_instance().get_numImages())
{

}

void MyModel::fromPrior()
{
	mag0 = -100. + 200.*randomU();
	for(size_t i=0; i<delta_mag.size(); i++)
		delta_mag[i] = tan(M_PI*(randomU() - 0.5));
}

double MyModel::perturb()
{
	int which = randInt(2);
	if(which == 0)
	{
		mag0 += 200.*pow(10., 1.5 - 6.*randomU())*randn();
		mag0 = mod(mag0 + 100., 200.) - 100.;
	}
	else
	{
		int which2 = randInt(delta_mag.size());
		double u = 0.5 + atan(delta_mag[which2])/M_PI;
		u += pow(10., 1.5 - 6.*randomU())*randn();
		u = mod(u, 1.);
		delta_mag[which2] = tan(M_PI*(u - 0.5));
	}

	return 0.;
}

double MyModel::logLikelihood() const
{
	return 0.;
}

void MyModel::print(std::ostream& out) const
{
	out<<mag0<<' ';
	for(size_t i=0; i<delta_mag.size(); i++)
		out<<delta_mag[i]<<' ';
}

string MyModel::description() const
{
	return string("");
}

