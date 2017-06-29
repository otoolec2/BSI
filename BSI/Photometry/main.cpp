/*

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
*/


#include <iostream>
#include "DNest4/code/DNest4.h" // This will work assuming you created the environment variable DNEST4_PATH
#include "MyModel.h"
#include "Data.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	
	Data::get_instance().load("Data/test_metadata.txt", "Data/test_image.txt");
	
	Sampler<MyModel> sampler = setup<MyModel>(argc, argv);
	sampler.run();
	
	return 0;
}

