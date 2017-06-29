import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from matplotlib.patches import Ellipse
import sys

startTime = datetime.now()

sample = np.loadtxt('Output/Output_Sample_Data.txt', usecols = ((0,)), skiprows = 3)
xc = np.loadtxt('Output/Output_Sample_Data.txt', usecols = ((2,)), skiprows = 3)
yc = np.loadtxt('Output/Output_Sample_Data.txt', usecols = ((3,)), skiprows = 3)
f = np.loadtxt('Output/Output_Sample_Data.txt', usecols = ((4,)), skiprows = 3)
a = np.loadtxt('Output/Output_Sample_Data.txt', usecols = ((5,)), skiprows = 3)

samples = []
SAMPLE = []
X = []
Y = []
F = []
A = []

# samples is used to investigate the number of samples contining a source within this box. 
# If this is the same as the total effective sample size, we can be extremely confident a source is in this region
for i in range(len(xc)):
	if sample[i] not in samples:
		samples.append(sample[i])
	if xc[i] >= 642 and xc[i] <= 644 and yc[i] >= 629 and yc[i] <= 631:
		SAMPLE.append(i)
		X.append(xc[i])
		Y.append(yc[i])
		F.append(f[i])
		A.append(a[i])

#print len(samples)
#print len(X), len(xc)
		
X_avg = sum(X)/len(X)
Y_avg = sum(Y)/len(Y)
F_avg = sum(F)/len(F)
A_avg = sum(A)/len(A)

print X_avg, Y_avg
print " "
sys.stdout.flush()

X_stddev = 0.
Y_stddev = 0.

for i in range(len(X)):
	X_stddev += ((X[i] - X_avg)**2)/len(X)
	Y_stddev += ((Y[i] - Y_avg)**2)/len(Y)

X_stddev = X_stddev**0.5
Y_stddev = Y_stddev**0.5

print X_stddev, Y_stddev
print " "
sys.stdout.flush()

F_stddev = 0.
A_stddev = 0.
for i in range(len(X)):
	F_stddev += ((F[i] - F_avg)**2)/len(F)
	A_stddev += ((A[i] - A_avg)**2)/len(A)

F_stddev = F_stddev**0.5
A_stddev = A_stddev**0.5

print F_avg, F_stddev 
print " "
print A_avg, A_stddev
print " "
sys.stdout.flush()

# Using a circle for simplification, but in reality an ellipse would account for different std dev's of x and y.
circle = plt.Circle((X_avg, Y_avg), 3*((X_stddev+Y_stddev)/2), color = 'r', alpha = 0.5)
fig, ax = plt.subplots()
ax.add_artist(circle)
plt.plot(X, Y, 'bo')
plt.plot(X_avg, Y_avg, 'r*', ms =20)
plt.xlabel("X")
plt.ylabel("Y")
plt.show()

plt.plot(SAMPLE, F, 'bo')
plt.xlabel("Sample Number")
plt.ylabel("Flux")
plt.savefig("Source_Fluxes.png")
plt.show()

plt.plot(SAMPLE, A, 'bo')
plt.xlabel("Sample Number")
plt.ylabel("Alpha")
plt.savefig("Source_Alphas.png")
plt.show()

print datetime.now() - startTime
sys.stdout.flush()