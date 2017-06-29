/*

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
*/


#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include "Data.h"
#include <cmath>


using namespace std;
using namespace DNest4;

MyModel::MyModel()
:objects(4, Data::get_instance().get_NMAX(), false, 
	MyConditionalPrior( Data::get_instance().get_x_min(), Data::get_instance().get_x_max(),
						Data::get_instance().get_y_min(), Data::get_instance().get_y_max(), 
						Data::get_instance().get_f_min(), Data::get_instance().get_f_max(), 
						Data::get_instance().get_a_min(), Data::get_instance().get_a_max()), 
	PriorType::uniform)
,image(Data::get_instance().get_ni(),
	vector<long double>(Data::get_instance().get_nj()))
{

}

void MyModel::from_prior(RNG& rng)
{
	objects.from_prior(rng);
	calculate_image();
	mu = Data::get_instance().get_mu(); 
	sigma = Data::get_instance().get_sigma(); 
}

void MyModel::calculate_image()
{
	// Get coordinate stuff from data
	const vector< vector<double> >& x = Data::get_instance().get_x_rays();
	const vector< vector<double> >& y = Data::get_instance().get_y_rays();
	
	double a = Data::get_instance().get_a();
	double b = Data::get_instance().get_b();
	double c = Data::get_instance().get_c();
	double d = Data::get_instance().get_d();
	
	// Diff
	bool update = objects.get_removed().size() == 0;

	// Components
	const vector< vector<double> >& components = (update)?(objects.get_added())
                                                 :(objects.get_components());
	if(!update)
	{
		// Zero the image
		image.assign(Data::get_instance().get_ni(),
			vector<long double>(Data::get_instance().get_nj(), 0.));
	}

	double xc, yc, f, alpha, beta;
	double xx, yy;
	
	beta = 2.5;
	
	for(size_t k=0; k<components.size(); ++k)		//Need to edit this to calculate the image using a Moffat25 PSF.
	{
		xc = components[k][0]; yc = components[k][1];	//My components vector should only have 4 components (poor word choice, I know), {xc, yc, f, alpha}
		f = components[k][2]; alpha = components[k][3];		//Recall my Mock Image calculation is Sum[f*Moffat25(x-xc, y-yc, alpha), {f, x, y}]
														
		for(size_t i=0; i<image.size(); i++)
		{
			for(size_t j=0; j<image[i].size(); j++)
			{
				xx = (x[i][j] - xc);
				yy = (y[i][j] - yc);
				
				image[i][j] += f * pow((1.*a + b*(xx*xx + yy*yy)/(c*alpha*alpha)) , -(beta*d));
			}
		}
	}
	
}

double MyModel::perturb(RNG& rng)
{
	double logH = 0.;
	
	if( rng.rand() <= 0.75 )		
	{
		logH += objects.perturb(rng);
		calculate_image();
	}
	else
	{
		mu = rng.rand_int(6)*Data::get_instance().get_mu(); 
		sigma = rng.rand_int(6)*Data::get_instance().get_sigma(); 
	}
	
	return logH;
}

double MyModel::log_likelihood() const
{
	const vector< vector<double> >& data = Data::get_instance().get_image();
	//const vector< vector<double> >& sig = Data::get_instance().get_sigma();

	double logL = 0.;
	double var;
	var = sigma*sigma;
	for(size_t i=0; i<data.size(); i++)
	{
		for(size_t j=0; j<data[i].size(); j++)
		{
			//var = sigma*sigma;// + sig[i][j]*sig[i][j];
			logL += -0.5*log(2.*M_PI*var) -0.5*pow(data[i][j] - (image[i][j] + mu), 2)/var;
		}
	}

	return logL;
	
}

void MyModel::print(std::ostream& out) const
{
	out<<setprecision(6);
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			out<<image[i][j]<<' ';
	out<<setprecision(10);
	out<<mu<<' '<<sigma<<' ';
	//out<<sigma<<' ';
	objects.print(out); out<<' ';		//fairly certain this is the posterior_sample.txt output
}

string MyModel::description() const
{
	return string("objects ");
}

