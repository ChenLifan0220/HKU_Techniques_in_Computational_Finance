#Assignment 2 
#Question 3.2
#Implied Volatility Calculation
#Chen Lifan_3035419171
#encoding:utf-8

'''
European call/put option
S - Stock price
K - Strike price
T - Time to maturity
t - Any time less than T
tau - Volatility
r - Risk free interest
q - repo rate
sigma_true - The true value of sigma
'''

import scipy.stats
from scipy.stats import norm
import scipy
import sympy as sympy
import numpy as np
import matplotlib.pyplot as plt
from math import *
import csv as csv
import pandas as pd
import time

def BlackScholes(CallPutFlag, S, K, t, T, tau, r, q):
    d1 = (log(float(S)/K)+(r-q)*(T-t))/(tau*sqrt(T-t))+(1/2)*tau*sqrt(T-t)
    d2 = d1-tau*sqrt(T-t)
    if CallPutFlag =='C':
        N1 = norm.cdf(np.float64(d1))
        N2 = norm.cdf(np.float64(d2))
        return S*exp(-q*(T-t))*N1-K*exp(-r*(T-t))*N2
    elif CallPutFlag == 'P':
        N3 = norm.cdf(np.float64(-d1))
        N4 = norm.cdf(np.float64(-d2))
        return K*exp(-r*(T-t))*N4-S*exp(-q*(T-t))*N3

def Vega(S,K,T,t,sigma,r,q):
	tau = sympy.Symbol('tau')
	y = S*sympy.exp(-q*(T-t))*sympy.sqrt(T-t)*1/sympy.sqrt(2*pi)*sympy.exp(-((sympy.log(S/K)+(r-q)*(T-t))/(tau*sympy.sqrt(T-t))+(1/2)*tau*sympy.sqrt(T-t))**2/2)
	yPrime = y.evalf(subs = {tau:sigma})
	return yPrime

def Implied_Volatility_Calculation(CallPutFlag,S,K,T,t,r,q,market_value):
    sigmahat = sqrt(2*abs((log(S/K)+(r-q)*(T-t))/(T-t)))
	#print('sigmahat',sigmahat)
    tol = 1e-5
    sigma = sigmahat
    if CallPutFlag == 'C':
        C_true = market_value
    elif CallPutFlag == 'P':
        P_true = market_value
    sigmadiff = 1
    n = 1
    nmax = 100
    if S-K*exp(-r*T)>0:
        CLower = S-K*exp(-r*T)
    else:
        CLower = 0

    if K*exp(-r*T)-S>0:
        PLower = K*exp(-r*T)-S
    else:
        PLower = 0

    while(sigmadiff >= tol and n < nmax):
        if CallPutFlag == 'C':
            C = BlackScholes('C',S,K,t,T,sigma,r,q)
            #print('C',C)
            if (C<S*exp(-q*(T-t))-K*exp(-r*(T-t)) or C>S*exp(-q*(T-t))):
                return 'NaN'
                break
            else:
                Cvega = Vega(S,K,T,t,sigma,r,q)
                #print('Cvega',Cvega)
                if round(Cvega,6) == 0.0:
                    return 'NaN'
                    break
                increment = (C - C_true)/Cvega
                #print('increment',increment)
        elif CallPutFlag == 'P':
            P = BlackScholes('P',S,K,t,T,sigma,r,q)
            #print('P',P)
            if (P>K*exp(-r*(T-t)) or P<(K-S*exp(-r*(T-t)))):
                return "NaN"
                break
            else:
                Pvega = Vega(S,K,T,t,sigma,r,q)
                #print('Pvega',Pvega)
                if round(Pvega,6) == 0.0:
                    return 'NaN'
                    break
                increment = (P - P_true)/Pvega
                #print('increment',increment)
        sigma = sigma - increment
        #print('sigma',sigma)
        n = n+1
        sigmadiff=abs(increment)
        #print('sigmadiff',sigmadiff)
    return sigma

def main():
    T = 8/365.0
    t = 0
    r = 0.04
    q = 0.2
    #readfile
    marketdatafile = open('marketdata.csv', 'r')
    marketdata_reader = csv.reader(marketdatafile)

    totalbid31 = 0
    totalask31 = 0
    no31 = 0
    totalbid32 = 0
    totalask32 = 0
    no32 = 0
    totalbid33 = 0
    totalask33 = 0
    no33 = 0
    for row in marketdata_reader:
        if row[1] == '510050':
            if row[0] <'2016-Feb-16 09:31:00':
                no31 = no31 + 1
                totalbid31 = totalbid31 + float(row[3])
                totalask31 = totalask31 + float(row[5])
            elif row[0] <'2016-Feb-16 09:32:00':
                no32 = no32 + 1
                totalbid32 = totalbid32 + float(row[3])
                totalask32 = totalask32 + float(row[5])
            elif row[0] <'2016-Feb-16 09:33:00':
                no33 = no33 + 1
                totalbid33 = totalbid33 + float(row[3])
                totalask33 = totalask33 + float(row[5])
    averagebid31 = totalbid31 / no31
    averageask31 = totalask31 / no31
    averagebid32 = totalbid32 / no32
    averageask32 = totalask32 / no32
    averagebid33 = totalbid33 / no33
    averageask33 = totalask33 / no33

    #print('averagebid31',averagebid31)
    #print('averageask31',averageask31)
    #print('averagebid32',averagebid32)
    #print('averageask32',averageask32)
    #print('averagebid33',averagebid33)
    #print('averageask33',averageask33)
    #Implied_Volatility_Calculation(CallPutFlag,S,K,T,t,r,q,market_value)
    #print(Implied_Volatility_Calculation('P',averageask31,2.15,8/365.0,t,r,q,0.2085))
    #print(Implied_Volatility_Calculation('P',averagebid31,2.15,8/365.0,t,r,q,0.1788))

    marketdata = pd.read_csv('marketdata.csv')
    instruments = pd.read_csv('instruments.csv')
    data = pd.merge(marketdata, instruments, on = ['Symbol'], how = 'left')
    #条件筛选取值
    data31 = data[(data['LocalTime']<'2016-Feb-16 09:31:00') & (data['Symbol']!= 510050)]
    data32 = data[(data['LocalTime']<'2016-Feb-16 09:32:00') & (data['LocalTime']>='2016-Feb-16 09:31:00') & (data['Symbol']!= 510050)]
    data33 = data[(data['LocalTime']<'2016-Feb-16 09:33:00') & (data['LocalTime']>='2016-Feb-16 09:32:00') & (data['Symbol']!= 510050)]

    #for 31.csv
    BidCframe31 = pd.DataFrame(columns = ['Strike','BidVolC'])
    BidPframe31 = pd.DataFrame(columns = ['Strike','BidVolP'])
    AskCframe31 = pd.DataFrame(columns = ['Strike','AskVolC'])
    AskPframe31 = pd.DataFrame(columns = ['Strike','AskVolP'])
    for indexs in data31.index:
        #根据索引，从每一行中取值
        bid_market_value = float(data31.loc[indexs,['Bid1']].values)
        #print('bid_market_value',bid_market_value)
        Strike = float(data31.loc[indexs,['Strike']].values)
        #print('Strike',Strike)
        Option_type = data31.loc[indexs,['OptionType']].values
        #print('optiontype',Option_type)
        ask_market_value = float(data31.loc[indexs,['Ask1']].values)
        if Option_type == 'C':
            BidVolC = Implied_Volatility_Calculation(Option_type,averagebid31,Strike,T,t,r,q,bid_market_value)
            #print(BidVol)
            AskVolC = Implied_Volatility_Calculation(Option_type,averageask31,Strike,T,t,r,q,ask_market_value)
            #print(AskVol)
            BidCframe31.loc[indexs] = {'Strike':Strike, 'BidVolC':BidVolC}
            AskCframe31.loc[indexs] = {'Strike':Strike, 'AskVolC':AskVolC}  
        elif Option_type == 'P':
            BidVolP = Implied_Volatility_Calculation(Option_type,averagebid31,Strike,T,t,r,q,bid_market_value)
            #print(BidVol)
            AskVolP = Implied_Volatility_Calculation(Option_type,averageask31,Strike,T,t,r,q,ask_market_value)
            #print(AskVol)
            BidPframe31.loc[indexs] = {'Strike':Strike, 'BidVolP':BidVolP}
            AskPframe31.loc[indexs] = {'Strike':Strike, 'AskVolP':AskVolP}
    #去重
    frame311 = BidCframe31.drop_duplicates(subset = ['Strike'], keep = 'last')
    #print(frame311)
    frame312 = AskCframe31.drop_duplicates(subset = ['Strike'], keep = 'last')
    #print(frame312)
    frame313 = BidPframe31.drop_duplicates(subset = ['Strike'], keep = 'last')
    #print(frame313)
    frame314 = AskPframe31.drop_duplicates(subset = ['Strike'], keep = 'last')
    #print(frame314)
    frame315 = pd.merge(frame313,frame314, on = ['Strike'], how = 'left')
    frame316 = pd.merge(frame311,frame312, on = ['Strike'], how = 'left')
    frame31 = pd.merge(frame315, frame316, on = ['Strike'], how = 'left')
    #print(frame31)
    frame_31 = frame31.sort_values(by = ['Strike'])
    filename = '31.csv'
    csvfile31 = open(filename,'w')
    writer = csv.writer
    frame_31.to_csv(filename, index = False, sep = ',')
    df31 = frame_31.set_index('Strike')
    #print(df31)
    plt.figure()
    plt.plot(df31['BidVolP'],label = 'BidVolP')
    plt.plot(df31['AskVolP'],label = 'AskVolP')
    plt.plot(df31['BidVolC'],label = 'BidVolC')
    plt.plot(df31['AskVolC'],label = 'AskVolC')
    plt.legend(loc = 'upper left') 
    plt.xlabel('Strike')
    plt.ylabel('Implied Volatility')
    plt.title('Plot_31')
    plt.savefig('31.png')
    plt.close()

    #for 32.csv
    BidCframe32 = pd.DataFrame(columns = ['Strike','BidVolC'])
    BidPframe32 = pd.DataFrame(columns = ['Strike','BidVolP'])
    AskCframe32 = pd.DataFrame(columns = ['Strike','AskVolC'])
    AskPframe32 = pd.DataFrame(columns = ['Strike','AskVolP'])
    for indexs in data32.index:
        bid_market_value = float(data32.loc[indexs,['Bid1']].values)
        Strike = float(data32.loc[indexs,['Strike']].values)
        Option_type = data32.loc[indexs,['OptionType']].values
        ask_market_value = float(data32.loc[indexs,['Ask1']].values)
        if Option_type == 'C':
            BidVolC = Implied_Volatility_Calculation(Option_type,averagebid32,Strike,T,t,r,q,bid_market_value)
            AskVolC = Implied_Volatility_Calculation(Option_type,averageask32,Strike,T,t,r,q,ask_market_value)
            BidCframe32.loc[indexs] = {'Strike':Strike, 'BidVolC':BidVolC}
            AskCframe32.loc[indexs] = {'Strike':Strike, 'AskVolC':AskVolC}  
        elif Option_type == 'P':
            BidVolP = Implied_Volatility_Calculation(Option_type,averagebid32,Strike,T,t,r,q,bid_market_value)
            AskVolP = Implied_Volatility_Calculation(Option_type,averageask32,Strike,T,t,r,q,ask_market_value)
            BidPframe32.loc[indexs] = {'Strike':Strike, 'BidVolP':BidVolP}
            AskPframe32.loc[indexs] = {'Strike':Strike, 'AskVolP':AskVolP}
    frame321 = BidCframe32.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame322 = AskCframe32.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame323 = BidPframe32.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame324 = AskPframe32.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame325 = pd.merge(frame323,frame324, on = ['Strike'], how = 'left')
    frame326 = pd.merge(frame321,frame322, on = ['Strike'], how = 'left')
    frame32 = pd.merge(frame325, frame326, on = ['Strike'], how = 'left')
    frame_32 = frame32.sort_values(by = ['Strike'])
    #print(frame32)
    filename = '32.csv'
    csvfile32 = open(filename,'w')
    writer = csv.writer
    frame_32.to_csv(filename, index = False, sep = ',')
    df32 = frame_32.set_index('Strike')
    #print(df32)
    plt.figure()
    plt.plot(df32['BidVolP'],label = 'BidVolP')
    plt.plot(df32['AskVolP'],label = 'AskVolP')
    plt.plot(df32['BidVolC'],label = 'BidVolC')
    plt.plot(df32['AskVolC'],label = 'AskVolC')
    plt.legend(loc = 'upper left') 
    plt.xlabel('Strike')
    plt.ylabel('Implied Volatility')
    plt.title('Plot_32')
    plt.savefig('32.png')
    plt.close()

    #for 33.csv
    BidCframe33 = pd.DataFrame(columns = ['Strike','BidVolC'])
    BidPframe33 = pd.DataFrame(columns = ['Strike','BidVolP'])
    AskCframe33 = pd.DataFrame(columns = ['Strike','AskVolC'])
    AskPframe33 = pd.DataFrame(columns = ['Strike','AskVolP'])
    for indexs in data33.index:
        bid_market_value = float(data33.loc[indexs,['Bid1']].values)
        Strike = float(data33.loc[indexs,['Strike']].values)
        Option_type = data33.loc[indexs,['OptionType']].values
        ask_market_value = float(data33.loc[indexs,['Ask1']].values)
        if Option_type == 'C':
            BidVolC = Implied_Volatility_Calculation(Option_type,averagebid33,Strike,T,t,r,q,bid_market_value)
            AskVolC = Implied_Volatility_Calculation(Option_type,averageask33,Strike,T,t,r,q,ask_market_value)
            BidCframe33.loc[indexs] = {'Strike':Strike, 'BidVolC':BidVolC}
            AskCframe33.loc[indexs] = {'Strike':Strike, 'AskVolC':AskVolC}  
        elif Option_type == 'P':
            BidVolP = Implied_Volatility_Calculation(Option_type,averagebid33,Strike,T,t,r,q,bid_market_value)
            AskVolP = Implied_Volatility_Calculation(Option_type,averageask33,Strike,T,t,r,q,ask_market_value)
            BidPframe33.loc[indexs] = {'Strike':Strike, 'BidVolP':BidVolP}
            AskPframe33.loc[indexs] = {'Strike':Strike, 'AskVolP':AskVolP}
    frame331 = BidCframe33.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame332 = AskCframe33.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame333 = BidPframe33.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame334 = AskPframe33.drop_duplicates(subset = ['Strike'], keep = 'last')
    frame335 = pd.merge(frame333,frame334, on = ['Strike'], how = 'left')
    frame336 = pd.merge(frame331,frame332, on = ['Strike'], how = 'left')
    frame33 = pd.merge(frame335, frame336, on = ['Strike'], how = 'left')
    frame_33 = frame33.sort_values(by = ['Strike'])
    #print(frame33)
    filename = '33.csv'
    csvfile31 = open(filename,'w')
    writer = csv.writer
    frame_33.to_csv(filename, index = False, sep = ',')
    df33 = frame_33.set_index('Strike')
    #print(df33)
    plt.figure()
    plt.plot(df33['BidVolP'],label = 'BidVolP')
    plt.plot(df33['AskVolP'],label = 'AskVolP')
    plt.plot(df33['BidVolC'],label = 'BidVolC')
    plt.plot(df33['AskVolC'],label = 'AskVolC')
    plt.legend(loc = 'upper left') 
    plt.xlabel('Strike')
    plt.ylabel('Implied Volatility')
    plt.title('Plot_33')
    plt.savefig('33.png')
    plt.close()


if __name__ == "__main__":
    main()
    # import profile
    # profile.run("main()")