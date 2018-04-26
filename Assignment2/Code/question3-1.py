#Assignment 2 
#Question 3.1
#Implied Volatility Calculation
#Chen Lifan_3035419171

'''
European call/put option

S = 2 - Stock price
K = 2 - Strike price
T = 3 - Time to maturity
t = 0 - Any time less than T
tau = 3 - Volatility
r = 0.03 - Risk free interest
q - repo rate
sigma_true = 0.3 - The true value of sigma
'''

import scipy.stats
from scipy.stats import norm
from sympy import *
import sympy as sympy
import numpy as np
import matplotlib.pyplot as plt
from math import *
import math as math


def BlackScholes(CallPutFlag,S,K,t,T,tau,r,q):
	d1 = (log(float(S)/K)+(r-q)*(T-t))/(tau*sqrt(T-t))+(1/2)*tau*sqrt(T-t)
	d2 = d1-tau*sqrt(T-t)
	x=Symbol('x')
	N1 = (1/sympy.sqrt(2*math.pi))*integrate(sympy.exp(-x**2/2),(x,float('-inf'),d1)).evalf()
	N2 = (1/sympy.sqrt(2*math.pi))*integrate(sympy.exp(-x**2/2),(x,float('-inf'),d2)).evalf()
	N3 = (1/sympy.sqrt(2*math.pi))*integrate(sympy.exp(-x**2/2),(x,float('-inf'),-d1)).evalf()
	N4 = (1/sympy.sqrt(2*math.pi))*integrate(sympy.exp(-x**2/2),(x,float('-inf'),-d2)).evalf()
	if CallPutFlag =='Call':
		return S*sympy.exp(-q*(T-t))*N1-K*sympy.exp(-r*(T-t))*N2
	elif CallPutFlag == 'Put':
		return K*sympy.exp(-r*(T-t))*N4-S*sympy.exp(-q*(T-t))*N3

def Vega(S,K,T,t,sigma,r,q):
	tau = Symbol('tau')
	y = S*sympy.exp(-q*(T-t))*sympy.sqrt(T-t)*1/sympy.sqrt(2*math.pi)*sympy.exp(-((sympy.log(S/K)+(r-q)*(T-t))/(tau*sympy.sqrt(T-t))+(1/2)*tau*sympy.sqrt(T-t))**2/2)
	yPrime = y.evalf(subs = {tau:sigma})
	return yPrime

def Implied_Volatility_Calculation(CallPutFlag,S,K,T,t,r,q,market_value):
	sigmahat = sqrt(2*abs((log(S/K)+(r-q)*(T-t))/(T-t)))
	#print('sigmahat',sigmahat)
	tol = 1e-8
	sigma = sigmahat
	if CallPutFlag == 'Call':
		C_true = market_value
	elif CallPutFlag == 'Put':
		P_true = market_value
	sigmadiff = 1
	n = 1
	nmax = 100

	while(sigmadiff >= tol and n < nmax):
		if CallPutFlag == 'Call':
			C = BlackScholes('Call',S,K,t,T,sigma,r,q)
			#print('C',C)
			if (C<S*exp(-q*(T-t))-K*exp(-r*(T-t)) or C>S*exp(-q*(T-t))):
				return np.NaN
				break
			else:
				Cvega = Vega(S,K,T,t,sigma,r,q)
				#print('Cvega',Cvega)
				if round(Cvega,6) == 0.0:
					return np.NaN
					break
				increment = (C - C_true)/Cvega
				#print('increment',increment)
		elif CallPutFlag == 'Put':
			P = BlackScholes('Put',S,K,t,T,sigma,r,q)
			#print('P',P)
			if (P>K*exp(-r*(T-t)) or P<(K-S*exp(-r*(T-t)))):
				return np.NaN
				break
			else:
				Pvega = Vega(S,K,T,t,sigma,r,q)
				#print('Pvega',Pvega)
				if round(Pvega,6) == 0.0:
					return np.NaN
					break
				increment = (P - P_true)/Pvega
				#print('increment',increment)
		sigma = sigma - increment
		#print('sigma',sigma)
		n = n+1
		sigmadiff=abs(increment)
		#print('sigmadiff',sigmadiff)
	return sigma

#print(BlackScholes('Put',2,2,0,3.0,0.03,0.35,0))
#print(Implied_Volatility_Calculation('Put',1.9582499999999998,2.25,8/365.0,0,0.04,0.2,0.2508))
#print(Implied_Volatility_Calculation('Put',1.95825,2.5,8/365.0,0,0.04,0.2,0.4828))

CallPutFlag = input("Option Type:")
S = input("Spot Price:")
K = input("K:")
T = input("T:")
t = input("t:")
r = input("r:")
q = input("q:")
market_value = input("market_value:")
print(Implied_Volatility_Calculation(CallPutFlag,float(S),float(K),float(T),float(t),float(r),float(q),float(market_value)))