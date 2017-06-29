/*

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
*/

// This is an adaptation of the header file Data.h found in DNest4/code/Examples/RJObect_GalaxyField.
// It requires relatively little adaptation as consistency with the data format used in the galaxy field example
// has been pursued, in the interest of avoiding the necessity of editing the structure of the RJObect template.

#ifndef DNest4_Photometry_Data
#define DNest4_Photometry_Data

#include <vector>

class Data
{
	private:
		// Number of pixels
		int ni, nj;
		
		//Maximum number of sources believed to be in image
		int NMAX;

		// Coordinates of image boundaries
		double x_min, x_max, y_min, y_max;
		
		//Limits of flux
		double f_min, f_max;
		
		//Limits of alpha
		double a_min, a_max;

		//parameters for Moffat profile
		double a, b, c, d;
		
		//mean and std dev of noise profile
		double mu, sigma;
		
		// Pixel widths
		double dx, dy;

		// Coordinates of pixel centers
		std::vector< std::vector<double> > x_rays;
		std::vector< std::vector<double> > y_rays;

		// The pixels
		std::vector< std::vector<double> > image;

		void compute_ray_grid();

	public:
		Data();
		void load(const char* metadata_file, const char* image_file);

		// Getters
		int get_ni() const { return ni; }
		int get_nj() const { return nj; }
		int get_NMAX() { return NMAX; }
		double get_x_min() const { return x_min; }
		double get_x_max() const { return x_max; }
		double get_y_min() const { return y_min; }
		double get_y_max() const { return y_max; }
		double get_f_min() const { return f_min; }
		double get_f_max() const { return f_max; }
		double get_a_min() const { return a_min; }
		double get_a_max() const { return a_max; }
		double get_dx() const { return dx; }
		double get_dy() const { return dy; }
		double get_a() const { return a; }
		double get_b() const { return b; }
		double get_c() const { return c; }
		double get_d() const { return d; }
		double get_mu() const { return mu; }
		double get_sigma() const { return sigma; }
		const std::vector< std::vector<double> >& get_x_rays() const
			{ return x_rays; }
		const std::vector< std::vector<double> >& get_y_rays() const
			{ return y_rays; }
		const std::vector< std::vector<double> >& get_image() const
			{ return image; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif

