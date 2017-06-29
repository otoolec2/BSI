# -*- coding: utf-8 -*-
"""
Created on Tue May 02 15:13:15 2017

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017

"""

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from datetime import datetime

def PSF(x, y, s1, s2, w):           # Will move this to another library at a later date, but for now, it can live here
    a = np.exp( -(x**2 + y**2)/( 2*(s1**2) ))
    b = np.exp( -(x**2 + y**2)/( 2*(s2**2) ))
    c = w/( 2.*np.pi*(s1**2) )
    d = ( 1. - w )/( 2.*np.pi*(s2**2) )
    
    return a*c + b*d                # See Brewer et al. for this definition of the PSF.

def Moffat_PSF(x, y, alpha):    # beta was originally a parameter, but for now we'll use a fixed value of 2.5

    beta = 2.5
    
    a = (x**2 + y**2)/alpha**2
    #b = (beta - 1)/np.pi*(alpha**2)
    
    return ((1 + a)**(-beta))

def Broken_Power_Law(x, h1, h2, alpha1, alpha2):
    Z1 = (alpha1**(-1))*(h1**(-alpha1) - h2**(-alpha1)) + (alpha2**(-1))*(h2**(-alpha1))
    Z2 = Z1*(h2**(alpha1-alpha2))
    
    if x < h1:
        return 0
    elif x >= h1 and x <= h2:
        return (Z1**(-1))*(x**(-alpha1 - 1))
    else:
        return (Z2**(-1))*(x**(-alpha2 - 1))

def MockImage(fluxes, coords, x, y, s1, s2, w):   # This yields the image we would expect to see in the absence of noise
    M = 0                                                   # Somewhat messy having to pass in all params for PSF, will try to think 
                                                          # of ways to tidy it up...
                                                          
    for i in range(len(fluxes)):
        M += fluxes[i]*PSF(x-coords[0][i], y-coords[1][i], s1, s2, w)
    
    return M                                                # Will have to run this for every (x, y) in image, so must determine pixel size
    
def MockImage_Moffat(fluxes, coords, x, y, alpha):   # This yields the image we would expect to see in the absence of 
    M = 0                                                 # noise, using a Moffat profile
                                                                     
    for i in range(len(fluxes)):
        M += fluxes[i]*Moffat_PSF(x-coords[0][i], y-coords[1][i], alpha)
    
    return M 

def Noise():      #Calculate noise level. Note, will have to do this at each x, y, but wanted to avoid passing in too many arguments
    fwhm = 0.06    
    
    std_dev = fwhm/( 2 * np.sqrt( 2 * np.log(2) ) )
    
    Epsilon = np.random.normal( 0.04, std_dev)
    
    return Epsilon

def Image(noise, MI):   #This is what the "noisy" image will look like
    
    I = MI + noise
    
    return I

def Coord(Xmin, Xmax, Xran, Ymin, Ymax, Yran, n):                           #Function to generate x, y coords for n stars.
    xcoords = np.random.uniform(Xmin - 0.1*Xran, Xmax + 0.1*Xran, n)        # Generate n random x,y coordinates for the stellar
    ycoords = np.random.uniform(Ymin - 0.1*Yran, Ymax + 0.1*Yran, n)        # sample, assuming a uniform distribution across the sky.
                                                                            # +/- 0.1xran/0.1yran allows for stars just outside the image
    return [xcoords, ycoords]                                               # but whose PSF bled into the image

def Write_to_FITS( twod_list, filename ):
    
    twodlist = np.array(twod_list)

    hduarr = fits.PrimaryHDU(twodlist)
    hdulist = fits.HDUList([hduarr])
    
    hdulist.writeto(filename + '.fits', clobber = True) # Note: According to the official astropy documentation
                                                        # http://docs.astropy.org/en/stable/io/fits/api/files.html#astropy.io.fits.writeto
                                                        # the keyword "overwrite" should be used, not "clobber", which appears to be an outdated keyword. 
                                                        # It may that the version of astropy available with the Anaconda package has not yet been updated, 
                                                        # as using "overwrite" returned an"Unexpected keyword" error. Updating the version of astropy on 
                                                        # our end should solve this, and users with the newer version will experience no issues. 
                                                        # Parking this for now though.

def Write_Data(FileName, n, in_frame, coords, intensities, a, b, pixels ): #May add S/N at a later time
    
    datafile = open(FileName + '.dat', 'w')

    header = "Parameters of the sample image generated in RandImage.py. \n \n"
    key = "N is the total number of stars in the sample. \nN_v is the number of stars within the frame of the image. \n" 
    key2 = "Alpha is the Alpha parameter of the Moffat profile describing the PSF. \nBeta is the Beta parameter of the Moffat Profile (assumed fixed at 2.5). \n \n"
    
    datafile.write(header)
    datafile.write(key)
    datafile.write(key2)
    datafile.write("Frame Size = " + str(pixels) + "x" + str(pixels) + " pixels \n \n")     # For now just assume a square array is used. Can generalize later.
    
    datafile.write("N = " + str(n) + "\n")
    datafile.write("N_v = " + str(in_frame) + "\n")
    datafile.write("Alpha = " + str(a) + "\n")
    datafile.write("Beta = " + str(b) + "\n \n")
    
    datafile.write("        X               Y              Flux       \n")
    datafile.write("--------------------------------------------------\n")

    for i in range(n):
        datafile.write( " " + str(coords[0][i]) + "  " + str(coords[1][i]) + "  " + str(intensities[i]) + "  \n")
    datafile.close()             

def Plot_Sample(PlotName, coords, xmin, ymin):               # Basic for now, might be added to later                                     
                    
    plt.figure(num = None, figsize = (6, 6))
    plt.plot(coords[0], coords[1], 'ko')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.savefig(PlotName + ".png")                               

    
    
startTime = datetime.now()

N = 2	    # Total number of stars in the random sample.
            # Note: Like Brewer et al, we allow the possibility of stars being outside the image proper, 
            # but bleeding in from the edges

# Placeholder for flux distribution
fluxes = [ np.random.uniform(0.0, 4.0) for k in range(N) ] #np.random.uniform(1.0e-3, 1.0e3, N) """np.random.uniform(0.0, 4.0)"""

xran = 2. # Range of x values
xmin = -1.  # Simply using Brewer et al.'s numbers for simplicity
xmax = xmin + xran

yran = 2. # Range of y values
ymin = -1.
ymax = ymin + yran

pixel_size = xran/10   # Want a 100x100 pixel array for now, but this can be easily changed if necessary 

Final_Image = [[None]*int(xran/pixel_size) for k in range(int(xran/pixel_size))]
Mock = [[None]*int(xran/pixel_size) for k in range(int(xran/pixel_size))]
epsilon = [[None]*int(xran/pixel_size) for k in range(int(xran/pixel_size))]

pos = Coord(xmin, xmax, xran, ymin, ymax, yran, N)         

In_Image = 0        # In_Image counts the number of stars that are actually in the image (ie. in (xmin, xmax), (ymin, ymax))
for i in range(N):  # In_image is typically around 60-70 for N=100, similar to the value given in Brewer et al.
    if pos[0][i] > xmin and pos[0][i] < xmax and pos[1][i] > ymin and pos[1][i] < ymax:
        In_Image += 1

print In_Image

alpha = np.random.uniform(pixel_size, 5*pixel_size)

x = xmin
y = ymin
for i in range(int(xran/pixel_size)):
    for j in range(int(yran/pixel_size)):
        Mock[i][j] = MockImage_Moffat(fluxes, pos, x, y, alpha)     # Image in absence of noise
        epsilon[i][j] = Noise()
        Final_Image[i][j] = Image(epsilon[i][j], Mock[i][j])
        x += pixel_size
    
    y += pixel_size
    x = xmin

    
Plot_Sample("Stellar_Sample", pos, xmin, ymin)

Write_Data( "Sample_Data", N, In_Image, pos, fluxes, alpha, 2.5, int(xran/pixel_size) )

Write_to_FITS(Final_Image, 'Image_Sample')
Write_to_FITS(Mock, 'Mock_Sample')
Write_to_FITS(epsilon, 'Noise_Sample')

print datetime.now() - startTime