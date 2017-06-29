"""
Created on Sat May 13 14:22:49 2017

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
"""

# Pre-processing script which will create all data files necessary for the provided .cc code to run within DNest4.
# Prior to running this script, one must possess a FITs file of the image to be analysed, as well as a FITS file
# containing the PSF generated using http://tinytim.stsci.edu/cgi-bin/tinytimweb.cgi . 
# The outputs of this script are two data files which contain all parameters needed by DNest4. Most are generated 
# using the FITS files, but there are 8 parameters in particular which will require direct user interaction: 
# xmin, xmax, ymin, ymax, the mean and standard deviation of the noise profile, and the minimum and maximum possible 
# value of the alpha parameter of the Moffat25 profile used. 

# These user-input defined parameters need not be perfect, merely reasonable estimates. The mean and standard deviation
# of the noise profile can be estimated by examining a subsection of the larger image from which your test image is drawn.
# The alpha parameter of the Moffat25 profile describes how wide the PSF is. In this current work, reasonable values were 
# were found to be 0.5 and 2.5 for the minimum and maximum, respectively. The distribution from which values of alpha are 
# drawn in DNest 4 will be wider anyway, so there should be minimal risk that one has over- or under-estimated the minimum
# and maximum possible value of alpha. However, reasonable estimates should reduce the spread of the DNest4 samples around 
# the actual value.
# The values of xmin, xmax, ymin, ymax should be obvious from one's data file. Here, defaults of xmin = -1.0, xmax = 1.0 and
# the same for ymin, ymax are used, and they should suffice as DNest4 will simply return estimates of source positions between 
# -1.0 and 1.0. However, they can easily be changed if desired

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import scipy.optimize as opt
from datetime import datetime
import sys

startTime = datetime.now()

#Function describing a Moffat profile which will be fitted to the tinytim PSF imported
def Moffat_PSF((x,y), alpha, a, b, c, d):    

	beta = 2.5
	mid = (b*(x**2. + y**2.))/(c*(alpha**2.))
	
	mid2 = ((a*1. + mid)**(-beta*d))
	
	return np.array(mid2.ravel())
	
	
# Note: All print statements in code used here are followed by flushing the output. This is due this work being done on a Windows machine using the Cygwin environment. 
print ">>>>Beginning pre-processing."
print " "
sys.stdout.flush()


# Name of image FITS file and PSF FITS file
image_filename = "F814W_chip2_region.fits"
psf_filename = "WFPC2_F814W.fits"


# Location of output data files and PSF file 
data_path = "Photometry/Data/"
image_path = "Images/"
psf_path = "PSFs/"


# User input code

# Entering the limits of the x and y values the subject image runs over
print ">>>>Please enter the range of x values your image runs over (on one line, seperated by a space)."
print ">>>>Note: One can use simple values such as -1.0 1.0, but this may make it difficult to compare DNest4's predictions with initial guesses etc."
sys.stdout.flush()

flag = 0

while flag == 0:
	xmin, xmax = raw_input().split(' ')
	xmin = float(xmin)
	xmax = float(xmax)   
	
	if xmin > xmax:
		xmid = xmin
		xmid2 = xmax
		
		xmin = xmid2
		xmax = xmid
		
		option = 0
		print ">>>>You appear to have entered xmin before xmax. They have been switched. Do you wish to change this? y/n"
		sys.stdout.flush()
		while option == 0:
			cont = raw_input()
			
			if cont == 'y':
				print ">>>>Re-enter the range of x values."
				sys.stdout.flush()
				break
			elif cont == 'n':
				option = 1
				flag = 1
			else:
				print ">>>>Not a valid option. Enter y or n."
				sys.stdout.flush()
				
	elif xmin == xmax:
		print ">>>>You have entered identical values for xmin and xmax. Please re-enter these values."
		sys.stdout.flush()
		
	else:
		flag = 1

print " "
print ">>>>Please enter the range of y values your image runs over (on one line, seperated by a space)."
print ">>>>Note: One can use simple values such as -1.0 1.0, but this may make it difficult to compare DNest4's predictions with initial guesses etc."
sys.stdout.flush()

flag = 0

while flag == 0:
	ymin, ymax = raw_input().split(' ')
	ymin = float(ymin)
	ymax = float(ymax)   
	
	if ymin > ymax:
		ymid = ymin
		ymid2 = ymax
		
		ymin = ymid2
		ymax = ymid
		
		option = 0
		print ">>>>You appear to have entered ymin before ymax. They have been switched. Do you wish to change this? y/n"
		sys.stdout.flush()
		while option == 0:
			cont = raw_input()
			
			if cont == 'y':
				print ">>>>Re-enter the range of y values."
				sys.stdout.flush()
				break
			elif cont == 'n':
				option = 1
				flag = 1
			else:
				print ">>>>Not a valid option. Enter y or n."
				sys.stdout.flush()
				
	elif ymin == ymax:
		print ">>>>You have entered identical values for ymin and ymax. Please re-enter these values."
		sys.stdout.flush()
		
	else:
		flag = 1
	
	

# Entering the maximum number of sources believed to be in the subject image.
print " "
print ">>>>Please enter the maximum number of sources believed to be in the image."
print ">>>>The larger the number you enter, the longer DNest34 will take to complete, so it may not be worthwhile to choose a number you know is far above the number of sources in the image."
sys.stdout.flush()

flag = 0

while flag == 0:
	NMAX = int(raw_input())
	
	if NMAX < 0:
		print ">>>>There cannot be a negative number of sources in your image."
		sys.stdout.flush()
	elif NMAX == 0:
		print ">>>>Limiting the possible number of sources to 0 is not recommended, as DNest4 will produce no meaningful output."
		sys.stdout.flush()
	else:
		flag = 1



# Entering the estimated minimum and maximum values that the alpha parameter(s) could be.
print " "
print ">>>>Please enter the estimated minimum and maximum values of the parameter alpha to be used in the Moffat profile (on one line, seperated by a space):"
sys.stdout.flush()

flag = 0

while flag == 0:
	amin, amax = raw_input().split(' ')
	amin = float(amin)
	amax = float(amax)
	
	if amin <= 0. or amax <= 0.:
		print ">>>>Alpha cannot be negative or zero. Re-enter the estimated values."
		sys.stdout.flush()
	elif amin > 0. and amax > 0. and amin > amax:
		amid = amin
		amid2 = amax
		
		amin = amid2
		amax = amid
		
		option = 0
		print ">>>>You appear to have entered amin before amax. They have been switched. Do you wish to change this? y/n"
		sys.stdout.flush()
		while option == 0:
			cont = raw_input()
			
			if cont == 'y':
				print ">>>>Re-enter the estimated values."
				sys.stdout.flush()
				break
			elif cont == 'n':
				option = 1
				flag = 1
			else:
				print ">>>>Not a valid option. Enter y or n."
				sys.stdout.flush()
	else:
		flag = 1


# Entering the noise profiles parameters, mean and std.deviation
print " "		
print ">>>>Please enter the estimated mean and standard deviation of the noise profile (on one line, seperated by a space):"
sys.stdout.flush()

flag = 0

while flag == 0:
	mu, sigma = raw_input().split(' ')
	mu = float(mu)
	sigma = float(sigma)
	
	if mu <= 0. or sigma <= 0.:
		print ">>>>The mean and standard deviation cannot be negative. Re-enter the estimated values."
		sys.stdout.flush()
	else:
		flag = 1


# Entering the value which is assigned to bad pixels in your subject image.
print " "
print ">>>>Please enter the value assigned to bad data points in your image: (eg. -200)"
sys.stdout.flush()

bad_data = float(raw_input())




# All parameter entries are done now. The remaining code carries out the actual pre-processing necessary.

psf_hdu_list = fits.open(psf_path + psf_filename)

length = psf_hdu_list[0].data.shape[0]	# tinytim returns square pixel arrays

x0 = int(length/2)
y0 = int(length/2)

x = np.arange(-30, 30, 60./length)
y = np.arange(-30, 30, 60./length)
x,y = np.meshgrid(x, y)

psf_input = psf_hdu_list[0].data.ravel()	# .ravel() necessary as scipy.optimize requires 1D lists.
psf_input = np.array(psf_input)

psf_hdu_list.close()

init = [1.0, 1.0, 1.0, 1.0, 1.0]	# initial guess

# Fitting a Moffat profile to the PSF passed in.
popt, pcov = opt.curve_fit(Moffat_PSF, (x,y), psf_input, p0 = init)

# It should be noted that alpha is used as a fitting parameter here, but in DNest4 is used as a parameter to be varied in each sample.





image_hdu_list = fits.open(image_path + image_filename)
image = image_hdu_list[0].data    # May be necessary to change index to '1', depending on the FITS file used													 
data = open(data_path + "test_image.txt", 'w')

image_hdu_list.close()	
	
# Outputting the image data in a format compatible with DNest4
# Index value to prevent issue where image was flipped in the vertical direction.
print "Writing image data file."
sys.stdout.flush()
for i in range(len(image)):
    for j in range(len(image[i])):
		if image[-i-1][j] == bad_data:
			data.write(str(mu) + " ")
		else:
			data.write(str(image[-i-1][j]) + " ")
    data.write("\n")

data.close()
print "Done."
sys.stdout.flush()		




fmin = np.amin(image)	# Flux limits for DNest4. These will not be the actual limits of the uniform distributions used, which will be random,  
fmax = np.amax(image)	# merely values to define the distributions from which those limits will be drawn. 

if fmin < 0.0:
	fmin = 0.0

# Outputting parameters necessary for DNest4 algorithm
print "Writing metadata file."
sys.stdout.flush()

# Metadata file reads:
# xpix ypix xmin xmax ymin ymax fmin fmax amin amax a b c d mu sigma
# where xpix, ypix are the number of x pixels and number of y pixels, respectively. a, b, c, d are the parameters used to fit the Moffat profile to the tinytim PSF.
# The others should be self evident based on earlier parts of the code.
metadata = open(data_path + "test_metadata.txt", 'w')
metadata.write(str(len(image)) + " " + str(len(image[i])) + " " + str(NMAX) + " " + str(xmin) + " " + str(xmax) + " " + str(ymin) + " " + str(ymax)\
				+ " " + str(fmin) + " " + str(fmax) + " " + str(amin) + " " + str(amax) \
				+ " " + str(popt[1]) + " " + str(popt[2]) + " " + str(popt[3]) + " " + str(popt[4]) \
				+ " " + str(mu) + " " + str(sigma))
    
metadata.close()

print "Done."
sys.stdout.flush()

print datetime.now() - startTime
sys.stdout.flush()
