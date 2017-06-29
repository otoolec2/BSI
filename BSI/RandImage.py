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

def Moffat_PSF(x, y, alpha):    # Moffat 25 profile generator

    beta = 2.5
    
    a = (x**2 + y**2)/alpha**2
    #b = (beta - 1)/np.pi*(alpha**2)
    
    return ((1 + a)**(-beta))
    
def MockImage_Moffat(fluxes, coords, x, y, alpha):   # This yields the image we would expect to see in the absence of noise
    M = 0                                                 
	
    for i in range(len(fluxes)):
        M += fluxes[i]*Moffat_PSF(x-coords[0][i], y-coords[1][i], alpha)
    
    return M 

def Noise():      # Calculate noise level. Note, will have to do this at each x, y
    fwhm = 0.06    
    
    std_dev = fwhm/( 2 * np.sqrt( 2 * np.log(2) ) )
    
    Epsilon = np.random.normal( 0.04, std_dev)
    
    return Epsilon

def Image(noise, MI):   # This is what the "noisy" image will look like
    
    I = MI + noise
    
    return I

def Coord(Xmin, Xmax, Xran, Ymin, Ymax, Yran, n):                           # Function to generate x, y coords for n stars.
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

def Write_Data(FileName, n, in_frame, coords, intensities, a, b, pixels ): # Writes data file containing parameters for comparison with DNest4 output
    
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

def Plot_Sample(PlotName, coords, xmin, ymin):                                                   
                    
    plt.figure(num = None, figsize = (6, 6))
    plt.plot(coords[0], coords[1], 'ko')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.savefig(PlotName + ".png")                               

    
    
startTime = datetime.now()

N = 2	    # Total number of stars in the random sample.
            # Note: We allow the possibility of stars being outside the image proper, 
            # but bleeding in from the edges

# Simple uniform flux distribution
fluxes = [ np.random.uniform(0.0, 4.0) for k in range(N) ]

xran = 2. # Range of x values
xmin = -1.  
xmax = xmin + xran

yran = 2. # Range of y values
ymin = -1.
ymax = ymin + yran

pixel_size = xran/10   # Want a 10x10 pixel array

Final_Image = [[None]*int(xran/pixel_size) for k in range(int(xran/pixel_size))]
Mock = [[None]*int(xran/pixel_size) for k in range(int(xran/pixel_size))]
epsilon = [[None]*int(xran/pixel_size) for k in range(int(xran/pixel_size))]

pos = Coord(xmin, xmax, xran, ymin, ymax, yran, N)         

In_Image = 0        # In_Image counts the number of stars that are actually in the image (ie. in (xmin, xmax), (ymin, ymax))
for i in range(N):  
    if pos[0][i] > xmin and pos[0][i] < xmax and pos[1][i] > ymin and pos[1][i] < ymax:
        In_Image += 1

#print In_Image

# Simple uniform distribution for alpha
alpha = np.random.uniform(pixel_size, 5*pixel_size) 


# Generate image
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

   
# Output
Plot_Sample("Stellar_Sample", pos, xmin, ymin)

Write_Data( "Sample_Data", N, In_Image, pos, fluxes, alpha, 2.5, int(xran/pixel_size) )

Write_to_FITS(Final_Image, 'Images/Image_Sample')
Write_to_FITS(Mock, 'Images/Mock_Sample')
Write_to_FITS(epsilon, 'Images/Noise_Sample')

print datetime.now() - startTime