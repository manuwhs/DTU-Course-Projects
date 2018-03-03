
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import graph_lib as gr
import matplotlib.colors as ColCon

w = 10  # Width of the images
h = 6   # Height of the images

def get_SMA(time_series, L):
    """ Outputs the aritmetic mean of the time series using 
    a rectangular window of size L"""
    window = np.ones((L,1))
    sM = np.convolve(time_series.flatten(),window.flatten(), mode = "full")
#    print sM.shape
    sM = sM
    sM = sM/np.sum(window)    # Divide so that it is the actual mean.
    
    # The first window has to be special, coz otherwise it would be just 0s
#    for i in range(L):
#        sM[i] = sM[i]*L/(i+1)
    sM[:L] = (window * sM[L]).flatten()  # Set the first ones equal to the first fully obtained stimator
    
    
    sM = sM[:-L+1]    # Remove the last values since they hare convolved with 0's as well
    sM = sM.reshape ((sM.size,1))
    return sM

def get_WMA(time_series, L):
    """ Outputs the aritmetic mean of the time series using 
    a linearly descendent window of size L"""
    
    window = np.cumsum(np.ones((L,1))) / L   
    
    sM = np.convolve(time_series.flatten().flatten(),window.flatten(), mode = "full")
    
#    print sM.shape
    sM = sM
    sM = sM/np.sum(window)    # Divide so that it is the actual mean.
    
    sM[:L] = (np.ones((L,1)) * sM[L]).flatten()  # Set the first ones equal to the first fully obtained stimator
    sM = sM[:-L+1]    # Remove the last values since they hare convolved with 0's as well
    return sM

def get_EMA(time_series, L, alpha = -1):
    
    if (alpha == -1):
        alpha = 2.0/(L+1)
        
    """ Outputs the exponential mean of the time series using 
    a linearly descendent window of size L"""
    
    window = np.ones((L,1))
    factor = (1 - alpha)
    for i in range(L):
        window[i] *= factor
        factor *= (1 - alpha)
#        TODO  : Nada puto funciona !!!!!!!!!!
#    print np.asarray(time_series).shape, window.shape
    
    sM = np.convolve(time_series[:,0],window[:,0], mode = "full")
    
#    gr.plot_graph([],window,["fe","de","de"], 1)
#    print sM.shape
    sM = sM
    sM = sM/np.sum(window)    # Divide so that it is the actual mean.
    
    sM[:L] = (np.ones((L,1)) * sM[L]).flatten()  # Set the first ones equal to the first fully obtained stimator
    sM = sM[:-L+1]    # Remove the last values since they hare convolved with 0's as well
    return sM


def TrainedMean(time_series, L):
    """ First it trains the data so that the prediction is maximized"""
    
    ### Training phase, we obtained the MSQE of the filter for predicting the next value
    time_series = time_series.flatten()
    Ns = time_series.size
    
    Xtrain, Ytrain = fu.windowSample(time_series, L)
    

    window = np.linalg.pinv((Xtrain.T).dot(Xtrain))
    window = window.dot(Xtrain.T).dot(Ytrain)
    window = np.fliplr([window])[0]
    
#    gr.plot_graph([],window,["fe","de","de"], 1)
    
    sM = np.convolve(time_series.flatten(),window.flatten(), mode = "full")
    
#    print sM.shape
    sM = sM
    sM = sM/np.sum(window)    # Divide so that it is the actual mean.
    
    sM[:L] = (np.ones((L,1)) * sM[L]).flatten()  # Set the first ones equal to the first fully obtained stimator
    sM = sM[:-L+1]    # Remove the last values since they hare convolved with 0's as well
    return sM

def get_TrCrMr (time_series, alpha = -1):
    """ Triple Cruce de la Muerte. Busca que las exponenciales 4, 18 y 40 se crucen
    para ver una tendencia en el mercado despues de un tiempo lateral """
    L1 = 4
    L2 = 18
    L3 = 40
    
    eM1 = get_EMA(time_series, L1, alpha)
    eM2 = get_EMA(time_series, L2, alpha)
    eM3 = get_EMA(time_series, L3, alpha)
    
    return np.array([eM1,eM2,eM3]).T

def get_HMA (time_series, L):
    """ Hulls Moving Average !! L = 200 usually"""
    WMA1 = get_WMA(time_series, L/2) * 2
    WMA2 = get_WMA(time_series, L)
    
    HMA = get_WMA(WMA1 - WMA2, np.sqrt(L))
    
    return HMA
    
def get_MAg (time_series, L, alpha = -1):
    """ Generalized Moving Average from Hull"""
    EMA1 = get_EMA(time_series, L/2, alpha) * 2
    EMA2 = get_EMA(time_series, L, alpha)
    
    EMA = get_EMA(EMA1 - EMA2, np.sqrt(L), alpha)
    
    return EMA

def get_MACD (time_series, Ls = 12, Ll = 26, Lsmoth = 9, alpha = -1):
    """ MACD """
    
    eMlong = get_EMA(time_series, Ll, alpha)
    eMshort = get_EMA(time_series, Ls, alpha)
    
    MACD = get_SMA(eMshort - eMlong, Lsmoth)
    
    return MACD

def get_momentum (time_series, N = 1):
    
    momentum = np.array(time_series[N:] - time_series[:-N])
    zero_vec = np.zeros((N,1))  # Add zero vector
    momentum = np.concatenate((zero_vec,momentum), axis = 0)
    
    return momentum


def get_RSI (time_series, N = 1):
    # Relative strength index
    momentum = np.array(time_series[N:] - time_series[:-N])
    zero_vec = np.zeros((N,1))  # Add zero vector
    momentum = np.concatenate((zero_vec,momentum), axis = 0)
    
    return momentum


    
def Bollinger_Bands (time_series, mean):
    return 1
    