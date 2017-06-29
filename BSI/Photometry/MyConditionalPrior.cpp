/*

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
*/


#include "MyConditionalPrior.h"
#include "DNest4/code/DNest4.h"
#include <cmath>

using namespace DNest4;

MyConditionalPrior::MyConditionalPrior(double x_min, double x_max,
					double y_min, double y_max, 
					double fluxlim_min, double fluxlim_max,
					double alphalim_min, double alphalim_max)	//will make alphalim_max/min be dependent on pixel_size
:x_min(x_min)
,x_max(x_max)
,y_min(y_min)
,y_max(y_max)
,fluxlim_min(fluxlim_min)
,fluxlim_max(fluxlim_max)
,alphalim_min(alphalim_min)
,alphalim_max(alphalim_max)
{

}

void MyConditionalPrior::from_prior(RNG& rng)
{
	fmax = 5 * fluxlim_min + rng.rand() * 5 * ( fluxlim_max - fluxlim_min );
	fmin = fluxlim_min + rng.rand() * ( fluxlim_max - fluxlim_min );

	amax = 5 * alphalim_min + 5 * rng.rand() * ( alphalim_max - alphalim_min );
	amin = alphalim_min + rng.rand() * ( alphalim_max - alphalim_min );
}

double MyConditionalPrior::perturb_hyperparameters(RNG& rng)
{
	double logH = 0.;
	
	int which = rng.rand_int(4);
	
	if( which == 0 ){
		fmax = 5 * fluxlim_min + rng.rand() * 5 * ( fluxlim_max - fluxlim_min );
		//fmax = mod(fmax, 5*fluxlim_max);
	}
	else if( which == 1 ){
		fmin = fluxlim_min + rng.rand() * ( fluxlim_max - fluxlim_min );
		//fmin = mod(fmin, 5*fluxlim_min);
	}
	else if( which == 2 ){
		amax = 5 * alphalim_min + 5 * rng.rand() * ( alphalim_max - alphalim_min );
		//amax = mod(amax, 5*alphalim_max);
	}
	else{
		amin = alphalim_min + rng.rand() * ( alphalim_max - alphalim_min );
		//amin = mod(amin, 5*alphalim_min);
	}

	return logH;
}

double MyConditionalPrior::log_pdf(const std::vector<double>& vec) const
{
	if( vec[0] < x_min || vec[0] > x_max ||
		vec[1] < y_min || vec[1] > y_max ||
		vec[2] < fmin || vec[2] > fmax ||
		vec[3] < amin || vec[3] > amax )
	return -1e300;
	
	double logp = 0.;
	
	logp = -log( fmax - fmin ) - log( amax - amin );
	
	return logp;
}

void MyConditionalPrior::from_uniform(std::vector<double>& vec) const
{
	vec[0] = x_min + ( x_max - x_min ) * vec[0];
	vec[1] = y_min + ( y_max - y_min ) * vec[1];
	vec[2] = fmin + ( fmax - fmin ) * vec[2];
	vec[3] = amin + ( amax - amin ) * vec[3];
}

void MyConditionalPrior::to_uniform(std::vector<double>& vec) const
{
	vec[0] = ( vec[0] - x_min )/( x_max - x_min );
	vec[1] = ( vec[1] - y_min )/( y_max - y_min );
	vec[2] = ( vec[2] - fmin )/( fmax - fmin );
	vec[3] = ( vec[3] - amin )/( amax - amin );
}

void MyConditionalPrior::print(std::ostream& out) const
{
	out<<fmax<<' '<<fmin<<' '<<amax<<' '<<amin<<' ';
}

