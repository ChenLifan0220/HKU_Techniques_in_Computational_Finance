#Assignment 2 
#Question 2.2
#Correlated normal random variables
#Chen Lifan_3035419171

import scipy.stats as norm
import numpy as np
import matplotlib.pyplot as plt
from math import *

sampleNo = 200
mu = 0
sigma = 1
np.random.seed(0)
X = np.random.normal(mu, sigma, sampleNo)
Y = np.random.normal(mu, sigma, sampleNo)
Z = 0.5*X+sqrt(1-0.5*0.5)*Y

rau = np.cov(X, Z)/(np.std(X) * np.std(Z))
#plt.hist(X, 30, normed=True)
#plt.show()
print(rau[0,1], '\n')


