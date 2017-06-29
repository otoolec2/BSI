/*

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
*/


#ifndef DNest4_Template_MyConditionalPrior
#define DNest4_Template_MyConditionalPrior

#include "DNest4/code/DNest4.h"

class MyConditionalPrior:public DNest4::ConditionalPrior
{
	private:
		
		//limits on x, y coords, as well as pixel size, (used to generate value for alpha)
		double x_min, x_max, y_min, y_max;
		double fluxlim_min, fluxlim_max;
		double alphalim_min, alphalim_max;
		
		//upper and lower limits for flux from an individual star
		double fmax, fmin;
		
		//alpha parameter limits for the Moffat25 profile (we're assuming the beta parameter is 2.5)
		double amax, amin;
	
		double perturb_hyperparameters(DNest4::RNG& rng);
		
	public:
		MyConditionalPrior( double x_min, double x_max, double y_min, double y_max,
					double fluxlim_min, double fluxlim_max,
					double alphalim_min, double alphalim_max);

		void from_prior(DNest4::RNG& rng);

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
		static const int weight_parameter = 1;
};

#endif

