
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os as os
import matplotlib.colors as ColCon
from scipy import spatial
import datetime as dt
from sklearn import linear_model

#########################################################3
############### BASIC MATH ##############################
##########################################################

# These funcitons expect a price sequences:
#     has to be a np.array [Nsamples][Nsec]

def get_return(price_sequences):
    # Get the return of the price sequences
    Nsam, Nsec = price_sequences.shape
    R = (price_sequences[1:,:] - price_sequences[0:-1,:])/ price_sequences[0:-1,:]
    # Add zero vector so that the length of the output is the same
    # as of the input
    zero_vec = np.zeros((1,Nsec))  
    R = np.concatenate((zero_vec,R), axis = 0)

    return R
    
def get_cumReturn(price_sequences):
    # Get the cumulative return of the price sequences
    returns = get_return(price_sequences)
    cR = np.cumsum(returns, axis = 0)
    return cR
    
def get_SharpR(Returns, axis = 0, Rf = 0):
    # Computes the Sharp ratio of the given returns
    # Rf = Risk-free return

    E_Return = np.mean(Returns,axis)
    std_Return = np.std(Returns,axis)
    SR = (E_Return- Rf)/std_Return
    return SR

def get_SortinoR(Returns, axis = 0):
    # Computes the Sortino ratio of the given returns
    E_Return = np.mean(Returns,axis)
    Positive_Returns = Returns[np.where(Returns < 0)]
    std_Return = np.std(Positive_Returns,axis)
    SR = E_Return/std_Return
    return SR

def get_covMatrix(returns):
    # We need to transpose it to fit the numpy standard
    covRet = np.cov(returns.T)
    return covRet
    
def get_corrMatrix(returns):
    # We need to transpose it to fit the numpy standard
    covRet = np.corrcoef(returns.T)
    return covRet
    
def get_linearRef(X,Y):
    ## Calculates the parameters of the linear regression
    
    Nsam, Ndim = X.shape
#    X = np.concatenate((np.ones((Nsam,1)), X ), axis = 1)
    # Create linear regression object
    regr = linear_model.LinearRegression()
    
    # Train the model using the training sets
    regr.fit(X, Y)

#    coeffs = np.array([regr.intercept_, regr.coef_])[0]
    
    coeffs = np.append(regr.intercept_, regr.coef_)
#    coeffs = np.concatenate((regr.coef_,regr.intercept_), axis = 1)
    return coeffs
    

def obtain_equation_line(Rf, Epoint, STDpoint):
    # Given the Rf and a portfolio point
    # This function calculates the equation of the line
    # Where the portfolio should be.
    P1 = [STDpoint, Epoint]
    P0 = [0, Rf]   # Origin point
    
    slope = (P1[1] - P0[1])/(P1[0] - P0[0])
    bias = slope * P0[0] + P0[1]  
    
    param = [bias, slope]
    
    return param
    
def get_TurnOver(w1,w2):
    # Gets the turn over between two porfolios
    # That is the absolute value of allocation we have to do.
    to = np.subtract(w1,w2)
    to = np.abs(to)
    to = np.sum(to)
    
    return to