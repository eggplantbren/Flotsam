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
,n_qso(1000)
,tau(Data::get_instance().get_numImages())
,y_qso(1000)
{

}

void MyModel::fromPrior()
{
	mag0 = -100. + 200.*randomU();

	delta_mag[0] = 0.;
	for(size_t i=1; i<delta_mag.size(); i++)
		delta_mag[i] = tan(M_PI*(randomU() - 0.5));

	tau_qso = exp(log(1.) + log(1E6)*randomU());
	beta_qso = exp(log(1E-3) + log(1E6)*randomU());
	for(size_t i=0; i<n_qso.size(); i++)
		n_qso[i] = randn();

	tau[0] = 0.;
	for(size_t i=1; i<tau.size(); i++)
		tau[i] = 10.*tan(M_PI*(randomU() - 0.5));

	assemble();
}

double MyModel::perturb()
{
	double logH = 0.;
	int which = randInt(6);

	if(which == 0)
	{
		mag0 += 200.*pow(10., 1.5 - 6.*randomU())*randn();
		mag0 = mod(mag0 + 100., 200.) - 100.;
	}
	else if(which == 1)
	{
		int which2 = 1 + randInt(delta_mag.size() - 1);
		double u = 0.5 + atan(delta_mag[which2])/M_PI;
		u += pow(10., 1.5 - 6.*randomU())*randn();
		u = mod(u, 1.);
		delta_mag[which2] = tan(M_PI*(u - 0.5));
	}
	else if(which == 2)
	{
		tau_qso = log(tau_qso);
		tau_qso += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		tau_qso = mod(tau_qso - log(1.), log(1E6)) + log(1.);
		tau_qso = exp(tau_qso);
	}
	else if(which == 3)
	{
		beta_qso = log(beta_qso);
		beta_qso += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		beta_qso = mod(beta_qso - log(1E-3), log(1E6)) + log(1E-3);
		beta_qso = exp(beta_qso);
	}
	else if(which == 4)
	{
		int which2 = 1 + randInt(tau.size() - 1);
		double u = 0.5 + atan(tau[which2]/10.)/M_PI;
		u += pow(10., 1.5 - 6.*randomU())*randn();
		u = mod(u, 1.);
		tau[which2] = 10.*tan(M_PI*(u - 0.5));
	}
	else if(which == 5)
	{
		double chance = pow(10., 0.5 - 4.*randomU());
		double scale = pow(10., 1.5 - 6.*randomU());
		bool full = randomU() <= 0.3;
		for(size_t i=0; i<n_qso.size(); i++)
		{
			if(randomU() <= chance)
			{
				if(full)
					n_qso[i] = randn();
				else
				{
					logH -= -0.5*pow(n_qso[i], 2);
					n_qso[i] += scale*randn();
					logH += -0.5*pow(n_qso[i], 2);
				}
			}
		}


	}

	assemble();

	return logH;
}

void MyModel::assemble()
{
	double alpha = exp(-1./tau_qso);
	y_qso[0] = mag0 + beta_qso/sqrt(1. - alpha*alpha)*n_qso[0];
	for(size_t i=1; i<y_qso.size(); i++)
		y_qso[i] = mag0 + alpha*(y_qso[i-1] - mag0) + beta_qso*n_qso[i];
}

double MyModel::logLikelihood() const
{
	double logL = 0.;

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	out<<mag0<<' ';
	for(size_t i=0; i<delta_mag.size(); i++)
		out<<delta_mag[i]<<' ';
	out<<tau_qso<<' '<<beta_qso<<' ';

	for(size_t i=0; i<tau.size(); i++)
		out<<tau[i]<<' ';

	for(size_t i=0; i<y_qso.size(); i++)
		out<<y_qso[i]<<' ';
}

string MyModel::description() const
{
	return string("");
}

