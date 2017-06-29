/*

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
*/


// This is based heavily on the Model object in the Galaxy Field example.

#ifndef DNest4_Photometry_MyModel
#define DNest4_Photometry_MyModel

#define _USE_MATH_DEFINES

#include "DNest4/code/DNest4.h"
#include "MyConditionalPrior.h"
#include <ostream>

class MyModel
{
	private:
		DNest4::RJObject<MyConditionalPrior> objects;

		// The model image
		std::vector< std::vector<long double> > image;
		void calculate_image();

		// How many steps since image was computed from scratch
		int staleness;
		
		//Noise standard deviation and mean
		double mu, sigma;

	public:
		// Constructor only gives size of params
		MyModel();

		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif

