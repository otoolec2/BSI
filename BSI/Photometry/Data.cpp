/*

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
*/

// This is an adaptation of the source file Data.cpp found in DNest4/code/Examples/RJObect_GalaxyField.
// It requires relatively little adaptation as consistency with the data format used in the galaxy field example
// has been pursued, in the interest of avoiding the necessity of editing the structure of the RJObect template.

#include "Data.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* metadata_file, const char* image_file)
{
	/*
	* First, read in the metadata
	*/
	fstream fin(metadata_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<metadata_file<<"."<<endl;
	fin>>ni>>nj;
	fin>> NMAX;
	fin>>x_min>>x_max>>y_min>>y_max;
	fin>>f_min>>f_max>>a_min>>a_max;
	fin>>a>>b>>c>>d;	//a, b, c, d are the parameters used in the Moffat profile fitted to the tinytim PSF used
	fin>>mu>>sigma;
	fin.close();

	// Make sure maximum > minimum
	if(x_max <= x_min || y_max <= y_min)
		cerr<<"# ERROR: strange input in "<<metadata_file<<"."<<endl;

	// Compute pixel widths
	dx = (x_max - x_min)/nj;
	dy = (y_max - y_min)/ni;

	// Check that pixels are square
	if(abs(log(dx/dy)) >= 1E-3)
		cerr<<"# ERROR: pixels aren't square."<<endl;

	compute_ray_grid();

	/*
	* Now, load the image
	*/
	fin.open(image_file, ios::in);
	if(!fin)
		cerr<<"# ERROR: couldn't open file "<<image_file<<"."<<endl;
	image.assign(ni, vector<double>(nj));
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			fin>>image[i][j];
	fin.close();

}

void Data::compute_ray_grid()
{
	// Make vectors of the correct size
	x_rays.assign(ni, vector<double>(nj));
	y_rays.assign(ni, vector<double>(nj));

	// Distance between adjacent rays
	double L = dx;

	for(size_t i=0; i < x_rays.size(); i++)
	{
		for(size_t j=0; j < x_rays[i].size(); j++)
		{
			x_rays[i][j] = x_min + (j + 0.5)*L;
			y_rays[i][j] = y_max - (i + 0.5)*L;
		}
	}
}

