
"""
Created on Thu Jun 01 10:55:25 2017

@author: Conor O'Toole
		 University College Dublin,
		 Ireland
		 
@supervisor: Dr. Morgan Fraser
			 University College Dublin,
			 Ireland
			 
			 
		 Supported by University College Dublin Seed Funding Scheme 2017
"""

# Post-processing script which brings together all post-processing carried out by the standard DNest4 Python scripts 
# show_results.py, display.py and make_plots.py.
# Can be edited as needed. 

# Inputs for this script are the standard .txt. results files that are produced by the DNest4 algorithm, as well as
# the image and metadata files used as input in the DNest4 algorithm. Outputs include the standard statistical data 
# which are plotted by DNest4's standard post-processing, a histogram showing the number of sources detected in each 
# posterior sample, a data file containing the number of sources and parameters of each source, the integral under  
# the PSF of each source, as well as the posterior probability of each sample.

# The code also requires the user to input a first guess of the position of the source(s)

from pylab import *
import numpy as np
import dnest4 as dn4
import os
import copy
from datetime import datetime
import sys

startTime = datetime.now()

# Moffat profile needed to calculate integral under PSF
def Moffat_PSF(x, y, alpha, a, b, c, d):    # beta was originally a parameter, but for now we'll use a fixed value of 2.5

    beta = 2.5
    
    mid = (b*(x**2. + y**2.))/(c*(alpha**2.))
    
    return ((a*1. + mid)**(-beta*d))

# Altered the standard postprocess script from DNest4. This is intended to improve efficiency in obtatining certain parameter values.
# The code in this function is not commented as extensively as in the remainder of the script, given that it is essentially a minor adaptation 
# the DNest4 author's code.
def postprocess_v2(temperature=1., numResampleLogX=1, plot=True, loaded=[], \
			cut=0., save=True, zoom_in=True, compression_bias_min=1., verbose=True,\
			compression_scatter=0., moreSamples=1., compression_assert=None, single_precision=False):
	if len(loaded) == 0:
		levels_orig = np.atleast_2d(dn4.my_loadtxt("levels.txt"))
		sample_info = np.atleast_2d(dn4.my_loadtxt("sample_info.txt"))
	else:
		levels_orig, sample_info = loaded[0], loaded[1]

	# Remove regularisation from levels_orig if we asked for it
	if compression_assert is not None:
		levels_orig[1:,0] = -np.cumsum(compression_assert*np.ones(levels_orig.shape[0] - 1))

	cut = int(cut*sample_info.shape[0])
	sample_info = sample_info[cut:, :]

	if plot:
		plt.figure(1)
		plt.plot(sample_info[:,0], "k")
		plt.xlabel("Iteration")
		plt.ylabel("Level")

		plt.figure(2)
		plt.subplot(2,1,1)
		plt.plot(np.diff(levels_orig[:,0]), "k")
		plt.ylabel("Compression")
		plt.xlabel("Level")
		xlim = plt.gca().get_xlim()
		plt.axhline(-1., color='g')
		plt.axhline(-np.log(10.), color='g', linestyle="--")
		plt.ylim(ymax=0.05)

		plt.subplot(2,1,2)
		good = np.nonzero(levels_orig[:,4] > 0)[0]
		plt.plot(levels_orig[good,3]/levels_orig[good,4], "ko-")
		plt.xlim(xlim)
		plt.ylim([0., 1.])
		plt.xlabel("Level")
		plt.ylabel("MH Acceptance")

	# Convert to lists of tuples
	logl_levels = [(levels_orig[i,1], levels_orig[i, 2]) for i in range(0, levels_orig.shape[0])] # logl, tiebreaker
	logl_samples = [(sample_info[i, 1], sample_info[i, 2], i) for i in range(0, sample_info.shape[0])] # logl, tiebreaker, id
	logx_samples = np.zeros((sample_info.shape[0], numResampleLogX))
	logp_samples = np.zeros((sample_info.shape[0], numResampleLogX))
	logP_samples = np.zeros((sample_info.shape[0], numResampleLogX))
	P_samples = np.zeros((sample_info.shape[0], numResampleLogX))
	logz_estimates = np.zeros((numResampleLogX, 1))
	H_estimates = np.zeros((numResampleLogX, 1))

	# Find sandwiching level for each sample
	sandwich = sample_info[:,0].copy().astype('int')
	for i in range(0, sample_info.shape[0]):
		while sandwich[i] < levels_orig.shape[0]-1 and logl_samples[i] > logl_levels[sandwich[i] + 1]:
			sandwich[i] += 1

	for z in range(0, numResampleLogX):
		# Make a monte carlo perturbation of the level compressions
		levels = levels_orig.copy()
		compressions = -np.diff(levels[:,0])
		compressions *= compression_bias_min + (1. - compression_bias_min)*np.random.rand()
		compressions *= np.exp(compression_scatter*np.random.randn(compressions.size))
		levels[1:, 0] = -compressions
		levels[:, 0] = np.cumsum(levels[:,0])

		# For each level
		for i in range(0, levels.shape[0]):
			# Find the samples sandwiched by this level
			which = np.nonzero(sandwich == i)[0]
			logl_samples_thisLevel = [] # (logl, tieBreaker, ID)
			for j in range(0, len(which)):
				logl_samples_thisLevel.append(copy.deepcopy(logl_samples[which[j]]))
			logl_samples_thisLevel = sorted(logl_samples_thisLevel)
			N = len(logl_samples_thisLevel)

			# Generate intermediate logx values
			logx_max = levels[i, 0]
			if i == levels.shape[0]-1:
				logx_min = -1E300
			else:
				logx_min = levels[i+1, 0]
			Umin = np.exp(logx_min - logx_max)

			if N == 0 or numResampleLogX > 1:
				U = Umin + (1. - Umin)*np.random.rand(len(which))
			else:
				U = Umin + (1. - Umin)*np.linspace(1./(N+1), 1. - 1./(N+1), N)
			logx_samples_thisLevel = np.sort(logx_max + np.log(U))[::-1]

			for j in range(0, which.size):
				logx_samples[logl_samples_thisLevel[j][2]][z] = logx_samples_thisLevel[j]

				if j != which.size - 1:
					left = logx_samples_thisLevel[j+1]
				elif i == levels.shape[0]-1:
					left = -1E300
				else:
					left = levels[i+1][0]

				if j != 0:
					right = logx_samples_thisLevel[j-1]
				else:
					right = levels[i][0]

				logp_samples[logl_samples_thisLevel[j][2]][z] = np.log(0.5) + dn4.classic.logdiffexp(right, left)

		logl = sample_info[:,1]/temperature

		logp_samples[:,z] = logp_samples[:,z] - dn4.classic.logsumexp(logp_samples[:,z])
		logP_samples[:,z] = logp_samples[:,z] + logl
		logz_estimates[z] = dn4.classic.logsumexp(logP_samples[:,z])
		logP_samples[:,z] -= logz_estimates[z]
		P_samples[:,z] = np.exp(logP_samples[:,z])
		H_estimates[z] = -logz_estimates[z] + np.sum(P_samples[:,z]*logl)

		if plot:
			plt.figure(3)

			plt.subplot(2,1,1)
			plt.plot(logx_samples[:,z], sample_info[:,1], 'k.', label='Samples')
			plt.plot(levels[1:,0], levels[1:,1], 'g.', label='Levels')
			plt.legend(numpoints=1, loc='lower left')
			plt.ylabel('log(L)')
			plt.title(str(z+1) + "/" + str(numResampleLogX) + ", log(Z) = " + str(logz_estimates[z][0]))
			# Use all plotted logl values to set ylim
			combined_logl = np.hstack([sample_info[:,1], levels[1:, 1]])
			combined_logl = np.sort(combined_logl)
			lower = combined_logl[int(0.1*combined_logl.size)]
			upper = combined_logl[-1]
			diff = upper - lower
			lower -= 0.05*diff
			upper += 0.05*diff
			if zoom_in:
				plt.ylim([lower, upper])
			xlim = plt.gca().get_xlim()

		if plot:
			plt.subplot(2,1,2)
			plt.plot(logx_samples[:,z], P_samples[:,z], 'k.')
			plt.ylabel('Posterior Weights')
			plt.xlabel('log(X)')
			plt.xlim(xlim)

	P_samples = np.mean(P_samples, 1)
	P_samples = P_samples/np.sum(P_samples)
	logz_estimate = np.mean(logz_estimates)
	logz_error = np.std(logz_estimates)
	H_estimate = np.mean(H_estimates)
	H_error = np.std(H_estimates)
	ESS = np.exp(-np.sum(P_samples*np.log(P_samples+1E-300)))

	errorbar1 = ""
	errorbar2 = ""
	if numResampleLogX > 1:
		errorbar1 += " +- " + str(logz_error)
		errorbar2 += " +- " + str(H_error)

	if verbose:
		print("log(Z) = " + str(logz_estimate) + errorbar1)
		sys.stdout.flush()
		print("Information = " + str(H_estimate) + errorbar2 + " nats.")
		sys.stdout.flush()
		print("Effective sample size = " + str(ESS))
		sys.stdout.flush()

	# Resample to uniform weight
	N = int(moreSamples*ESS)
	w = P_samples
	w = w/np.max(w)
	rows = np.empty(N, dtype="int64")
	for i in range(0, N):
		while True:
			which = np.random.randint(sample_info.shape[0])
			if np.random.rand() <= w[which]:
				break
		rows[i] = which + cut

    # Get header row
	f = open("sample.txt", "r")
	line = f.readline()
	if line[0] == "#":
		header = line[1:]
	else:
		header = ""
	f.close()

	sample = dn4.loadtxt_rows("sample.txt", set(rows), single_precision)
	posterior_sample = None
	posterior_sample_2 = None
	if single_precision:
		posterior_sample = np.empty((N, sample["ncol"]), dtype="float32")
		posterior_sample_2 = np.empty((N, sample["ncol"] + 1), dtype="float32")
	else:
		posterior_sample = np.empty((N, sample["ncol"]))
		posterior_sample_2 = np.empty((N, sample["ncol"] + 1))

	for i in range(0, N):
		posterior_sample[i, :] = sample[rows[i]]	# This sets row i of posterior_sample to be that of sample.
		posterior_sample_2[i] = np.append(posterior_sample[i], sample_info[rows[i], 1])	# Adds the log likelihood to the end of the row
		

	if save:
		np.savetxt('weights.txt', w)
		if single_precision:
			np.savetxt("posterior_sample.txt", posterior_sample_2, fmt="%.7e",\
													header=header)
		else:
			np.savetxt("posterior_sample.txt", posterior_sample_2,\
													header=header)

	if plot:
		plt.show()

	return [logz_estimate, H_estimate, logx_samples]

# Piecewise linear stretch, necessary to resize image
def stretch(x):
	if( x.min() != x.max()):
		y = x.copy()
		y = (y - y.min())/(y.max() - y.min())
		y[y > 0.1] = 0.1 + 0.05*(y[y > 0.1] - 0.1)
	else:
		y = x.copy()
	return y
	


print ">>>>Beginning post-processing."
print " "
sys.stdout.flush()

output_header = "Output/"
data_header = "Photometry/Data/"

metadata_file = "Photometry/Data/test_metadata.txt"
	
# User input code

flag = 0
print " "
sys.stdout.flush()
print ">>>>Do you wish to view the statistical output of DNest4? y/n"
sys.stdout.flush()

while flag == 0:
	view_stats = raw_input()
	print " "
	sys.stdout.flush()
	
	if view_stats == 'y' or view_stats == 'n':
		flag = 1
		
	else:
		print ">>>>That is not a valid option. Please enter 'y' or 'n'."
		sys.stdout.flush()		

	
# Enter the number of sources believed to be in the image prior to running DNest4. 
# The estimated positions of each will then be requested, in order to plot them against the 
# predicted positions of sources.
print ">>>>Please enter the number of sources that is believed to be in the image."
print ">>>To skip this step, simply enter 0."
sys.stdout.flush()
flag = 0
while flag == 0:
	Ns = int(raw_input())
	
	if Ns < 0.:
		print "Image cannot contain a negative number of sources."
		sys.stdout.flush()
	else:
		flag = 1

if Ns > 0:
	print "Please enter the estimated positions of the {0} sources believed to be present in the image.".format(Ns)
	sys.stdout.flush()

	i = 1
	est_x = [[None] for k in range(Ns)]
	est_y = [[None] for k in range(Ns)]
	while i <= Ns:
		print "Source {0} estimated position (x and y coordinates seperated by a space).".format(i)
		sys.stdout.flush()
		est_x[i-1], est_y[i-1] = raw_input().split(' ')
		est_x[i-1] = float(est_x[i-1])
		est_y[i-1] = float(est_y[i-1])
		
		if est_x[i-1] < xmin or est_y[i-1] > xmax or est_y[i-1] < ymin or est_y[i-1] > ymax:
			print "Source {0} is outside the bounds of the image. Please re-enter the coordinates of this source.".format(i)
			sys.stdout.flush()
		else:
			i += 1

			
print " "
print ">>>>Do you wish to calculate the integrals under the PSF for each source in each sample? (y/n)"
print ">>>>WARNING: This will take a significant length of time. We recommend it only if hours of computer time can be spared."
sys.stdout.flush()

flag = 0
while flag == 0:
	integral_flag = raw_input()
	print " "
	sys.stdout.flush()
	
	if integral_flag == 'y' or integral_flag == 'n':
		flag = 1
		
	else:
		print ">>>>That is not a valid option. Please enter 'y' or 'n'."
		sys.stdout.flush()


print " "
print ">>>>Do you wish to generate images using the parameters in each sample? (y/n)"
print ">>>>WARNING: This will take a significant length of time. We recommend it only if hours of computer time can be spared."
print ">>>>If creating these images, do not forget to either move or delete images from a previous case."
sys.stdout.flush()
flag = 0
while flag == 0:
	picture_flag = raw_input()

	if picture_flag == 'y' or picture_flag == 'n':
		flag = 1
	
	else:
		print ">>>>That is not a valid option. Please enter 'y' or 'n'."
		sys.stdout.flush()






# Carry out the general post-processing from DNest4. This produces the effective sample file: posterior_sample.txt.
# These are all samples (usually about 10000) with a posterior weight above a certain threshold (ie. the most important samples)	

if view_stats == 'y':
		
		os.chdir('Photometry')
		Outputs = postprocess_v2(plot = True)
		os.chdir('../')
else:
	
	os.chdir('Photometry')
	Outputs = postprocess_v2(plot = False)							
	os.chdir('../')		
		
# Read in posterior_sample.txt and the image input used in DNest4. If one had a data file with the std deviation from the mean at each pixel,
# one can also read that in here. This was never necessary in our work, but perhaps will be for some. To include one, simply uncomment the 
# relevant line.
posterior_sample = atleast_2d(dn4.my_loadtxt('Photometry/posterior_sample.txt', single_precision=True))
sample = atleast_2d(dn4.my_loadtxt('Photometry/sample.txt'))
sample_info = atleast_2d(dn4.my_loadtxt('Photometry/sample_info.txt'))
data = loadtxt('Photometry/Data/test_image.txt')
#sig = loadtxt('Photometry/Data/test_sigma.txt')


# Read in the necessary values from the metadata file
# Conversion to integers necessary for some of the array operations which follow (DNest4 requires floats)
x_pix = np.loadtxt(metadata_file, usecols = (0,))
x_pix = int(x_pix)
y_pix = np.loadtxt(metadata_file, usecols = (1,))
y_pix = int(y_pix)
N_lim = np.loadtxt(metadata_file, usecols = (2,))
N_lim = int(N_lim)

xmin = np.loadtxt(metadata_file, usecols = (3,))
xmax = np.loadtxt(metadata_file, usecols = (4,))
ymin = np.loadtxt(metadata_file, usecols = (5,))
ymax = np.loadtxt(metadata_file, usecols = (6,))

# Moffat25 parameters necessary to calculate integral under PSF.
a = np.loadtxt(metadata_file, usecols = (11,))
b = np.loadtxt(metadata_file, usecols = (12,))
c = np.loadtxt(metadata_file, usecols = (13,))
d = np.loadtxt(metadata_file, usecols = (14,))




# Total number of pixels
N_sq = x_pix*y_pix

# Lists of the DNest4 parameters for each sample in the posterior sample.
# The number indicates the number of sources predicted by each sample, the fmax, fmin, amax, amin the limits of 
# the relevant uniform ditributions used in generating a given sample and the log_ls the log likelihood of that sample.
number = posterior_sample[:,N_sq + 8]
fmax = posterior_sample[:,N_sq + 4]
fmin = posterior_sample[:,N_sq + 5]
amax = posterior_sample[:,N_sq + 6]
amin = posterior_sample[:,N_sq + 7]
log_ls = posterior_sample[:,-1]

# The first entry in these following lists will contain the relevant parameters for the case where DNest 4 predicted zero sources.
# The subsequent entries will contain the data for each sample, and each source within those samples. These will not be grouped together
# based on number of sources in the way that the 0 source case is. This is because the parameters of a source in two sample which predict
# one source, say, can still differ dramatically. The overall probability of a given number of sources will be calculated later.
sample = [0]
num = [0.]
xc = [0.]
yc = [0.]
f = [0.]
alpha = [0.]
prob = [0.]

# In the following posterior probability calculation, we calculate the prior probability that each of the 4 parameters is within +- 0.0001 of the specific value.
# Thus, the posterior probability given is the prob that the parameters which would produce this data are within +- 0.0001 of these specific values.

# Calculate prior pdfs for x and y coordinates, as it will be the same for all samples.
xprior = 1./(xmax-xmin)
yprior = 1./(ymax-ymin)

print " "
print ">>>>Calculating posterior probabilities."
sys.stdout.flush()
for i in range(len(number)):
	prior = (0.0002**4)*np.abs(xprior) * np.abs(yprior) * (1./(np.abs(fmax[i] - fmin[i]))) * (1./(np.abs(amax[i] - amin[i]))) 
	posterior = prior * np.exp(log_ls[i] - Outputs[0])
	
	if number[i] > 0:
		for j in range(int(number[i])):
			sample.append(i+1)
			num.append(number[i])
			xc.append(posterior_sample[:,N_sq + 9 + j][i])
			yc.append(posterior_sample[:,N_sq + N_lim + 9 + j][i])
			f.append(posterior_sample[:,N_sq + 2*N_lim + 9 + j][i])
			alpha.append(posterior_sample[:,N_sq + 3*N_lim + 9 + j][i])
			prob.append(posterior)
	else:
		prob[0] += posterior 

### PENDING APPROVAL ###
P = 0
for i in range(len(num)):
	if num[i] == 0:
		P += prob[i]
	else:
		P += prob[i]/num[i]		# This accounts for multiple sources in the same sample

prob_temp = prob
for i in range(len(prob_temp)):
	prob_temp[i] = prob_temp[i]/P

prob = prob_temp

########################

print ">>>>Done."
print " "
sys.stdout.flush()

# Sum probabilities of samples with the same number of sources. The case with zero sources has already been done.
# The jth entry will  contain the total probability of j sources being present
NMAX = int(max(num))

Prob = [0. for k in range(int(NMAX) + 1)]
Prob[0] = prob[0] 

for i in range(1, len(sample)):
	
	for j in range(1, len(Prob)):
	
		if int(num[i]) == j:
		
			Prob[j] += prob[i]/float(j) # Division accounts for repeats from the same sample



			
# Calculating the integral under the PSF for each predicted source in each sample.
# This is a computationally expensive exercise, so if unnecessary can be omitted.

quarter_mark = int(len(sample)/4)

if integral_flag == 'y':
	
	print ">>>>Calculating integrals under PSFs."
	sys.stdout.flush()
	
	integral = []
	dx = 0.1
	dy = 0.1
	
	for i in range(len(sample)):
		
		if i == quarter_mark:
		
			print ">>>>25% done."
			sys.stdout.flush()
		elif i == 2*quarter_mark:
		
			print ">>>>50% done."
			sys.stdout.flush()
		elif i == 3*quarter_mark:
		
			print ">>>>75% done."
			sys.stdout.flush()
			
		
		integ = 0.
		X = xmin - np.abs(0.2 * xmin)	# abs accounts for cases where xmin is negative
		
		if num[i] != 0:
		
			while X <= xmax + 0.2 * xmax:
			
				Y = ymin - np.abs(0.2 * ymin)
				
				while Y <= ymax + 0.2 * ymax:
				
					integ += dx*dy*f[i]*Moffat_PSF(X-xc[i], Y - yc[i], alpha[i], a, b, c, d)
					Y += dy
				
				X += dx
		else:
		
			integ = 0.
			
		
		integral.append(integ)
		
		
	print ">>>>Done."





# Write the data for each sample. Samples containing multiple sources have that same number of lines.
# These lines share the same sample number, number of sources and posterior probability but each line contains the parameters for a
# different source in the sample.
print " "
print ">>>>Writing samples."
sys.stdout.flush()
sample_data = open(output_header + 'Output_Sample_Data.txt', 'w')
sample_data.write("#Sample_ID Number_of_Sources X Y Flux Alpha Flux_Integ Posterior \n \n")

if integral_flag =='y':

	for i in range(len(sample)):
		
		sample_data.write(str(sample[i]) + " " + str(num[i]) + " " + str(xc[i]) + " " + str(yc[i]) + " " + str(f[i]) + " " + str(alpha[i]) + " " + str(integral[i])+ " " + str(prob[i]) + "\n")

	sample_data.close()
else:

	for i in range(len(sample)):
		
		sample_data.write(str(sample[i]) + " " + str(num[i]) + " " + str(xc[i]) + " " + str(yc[i]) + " " + str(f[i]) + " " + str(alpha[i]) + " " + str(0.0) + " " + str(prob[i]) + "\n")

	sample_data.close()

print ">>>>Done"



# Write the metadata, including the effective sample size, the probabilities of a given number of sources being present.
print " "
print ">>>>Writing metadata."
sys.stdout.flush()

sample_metadata = open(output_header + 'Output_Sample_MetaData.txt', 'w')
sample_metadata.write("Effective Sample Size = " + str(len(number)) + "\n \n")
sample_metadata.write("#Number_of_sources Probability \n \n")

for j in range(len(Prob)):

	sample_metadata.write( str(j) + " " + str(Prob[j]) + "\n")

sample_metadata.close()
print ">>>>Done."




print " "
print ">>>>Producing plots."
sys.stdout.flush()

#Plot histogram of the number of sources predicted in each sample.
plt.hist(number, NMAX) 	
plt.xlabel('Number of Stars $N$', fontsize=14)
plt.ylabel('Number of posterior samples', fontsize=14)
plt.xlim([-0.5, number.max() + 1])
plt.savefig(output_header + 'N_star_result.png', bbox_inches='tight')
#plt.show()
plt.close()

# Plots of positions of sources from all samples with the same number of sources, versus the initial guesses entered by the user previously. 
# In addition, plots of fluxes and alpha parameters for samples with identical numbers of sources. Note, in the case where samples with more 
# than one source are considered (for instance, N=2), then for each point on the x-axis (sample number), there will be multiple points, one 
# for each source. No further processing is done here to calculate averages, etc, as DNest4 can mix up the order of sources between different 
# samples. Thus, this is left up to the user, if they determine from these plots that the spreads of the various parameters are narrow enough 
# to specify ranges in which the all of a specific parameter could be considered to have some from the same source. This will need to be determined
# on a case-by-case basis.
n = 1

while n <= NMAX:
	Xs = []
	Ys = []
	Fs = []
	As = []
	XAxis = []		# Array to plot on x axis, allowing plots of the flux and alpha values to be plotted
	for j in range(len(num)):
		if num[j] == n:	# The plots in cases where there are no sources are obvious.
			Xs.append(xc[j])
			Ys.append(yc[j])
			Fs.append(f[j])
			As.append(alpha[j])
			XAxis.append(j)
	
	if len(Xs) > 0:			# Only produce plots in those cases where at least one sample predicted n sources.
		plt.plot(Xs, Ys, 'bo')
		if Ns >0: 
			plt.plot(est_x, est_y, 'r*', ms = 15)
		plt.xlim((xmin, xmax))
		plt.ylim((ymin, ymax))
		plt.xlabel("X")
		plt.ylabel("Y")
		plt.savefig(output_header + 'FirstGuess_v_{0}_Source_Samples.png'.format(n))
		#plt.show()
		plt.close()
		
		plt.plot(XAxis, Fs, 'bo')
		plt.xlabel("Sample Number")
		plt.ylabel("Flux")
		plt.savefig(output_header + 'Flux_Distrib_{0}_Source Samples.png'.format(n))
		#plt.show()
		plt.close()
		
		plt.plot(XAxis, As, 'bo')
		plt.xlabel("Sample Number")
		plt.ylabel("Alpha")
		plt.savefig(output_header + 'Alpha_Distrib_{0}_Source Samples.png'.format(n))
		#plt.show()
		plt.close()
	
	n += 1

print ">>>>Done."




# The following plots a pixellated image with the predicted number of sources and their parameters for each sample. 
# This is a computationally expensive exercise, so if unnecessary can be omitted.
# If creating these images, do not forget to either move or delete images from a previous case, as if the new simulation has 
# less samples, images from previous instances will remain in the Frames file, potentially causing issues with analysis.

if picture_flag == 'y':

	print ">>>>Producing images based on the number of sources and their respective parameters for each sample."
	sys.stdout.flush()

	data = array(data)
	sig = array([[0.]*x_pix for k in range(y_pix)])
	#sig = array(sig)


	# For each sample in posterior_sample.txt, plot the image produced by the particular parameter values. 
	# This will act as essentially a catalogue of images approximate to the real thing
	# Will also plot the standardised residuals ((image-input_data)/sigma^2)

	for i in range(0, posterior_sample.shape[0]):
		
		img = posterior_sample[i, 0:N_sq].reshape((x_pix, y_pix)) # starts at 0, and uses the first x_pix*y_pix entries, where the image data is saved
		plt.subplot(1, 2, 1)
		plt.imshow(stretch(img), cmap='gray', interpolation='none')	# This plots the image generated from the posterior sample
		plt.title('Model {i}'.format(i=i+1))
		plt.gca().set_xticks([-0.0, (int(x_pix/2)+0.5), (x_pix - 0.5)])
		plt.gca().set_yticks([-0.0, (int(y_pix/2)+0.5), (y_pix - 0.5)])
		plt.gca().set_xticklabels(['-1', '0', '1'])
		plt.gca().set_yticklabels(['1', '0', '-1'])	# All images used in our work were assigned coordinates with limits [-1, 1], [1, 1]. If this is an issue, a simple relabelling will fix it
		plt.colorbar()
		
		plt.subplot(1, 2, 2)
		mu = posterior_sample[i, N_sq]
		sigma = sqrt(sig**2 + posterior_sample[i,N_sq + 1]**2)
		plt.imshow(-((img+mu) - data)/sigma, cmap='gray', interpolation='none')
		plt.title('Standardised Residuals')
		plt.gca().set_xticks([-0.0, (int(x_pix/2)+0.5), (x_pix - 0.5)])
		plt.gca().set_yticks([-0.0, (int(y_pix-2)+0.5), (y_pix - 0.5)])
		plt.gca().set_xticklabels(['-1', '0', '1'])
		plt.gca().set_yticklabels(['1', '0', '-1'])
		plt.colorbar()

		plt.savefig(output_header + 'Frames/' + '%0.4d'%(i+1) + '.png', bbox_inches='tight')
		#print('Output/Frames/' + '%0.4d'%(i+1) + '.png')

		#plt.show()		#Note, even if you are producing these plots, it is not recommended that you uncomment this line
		plt.close()
	print ">>>>Done."
	sys.stdout.flush()


print datetime.now() - startTime
sys.stdout.flush()