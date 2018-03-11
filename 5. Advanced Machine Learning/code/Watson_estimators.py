# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 20:47:19 2017

@author: montoya
"""

import import_folders
from scipy.special import hyp1f1
from scipy.special import gamma
from scipy.optimize import newton
import numpy as np
import utilities_lib as ul
import HMM_libfunc2 as HMMl

import Watson_distribution as Wad
import Watson_sampling as Was
import Watson_estimators as Wae
import general_func as gf
import warnings

def get_kappaNewton(k, args):  # The r is computed outsite
    Ndim = args[0]
    r = args[1]
    
    a = 0.5
    c = float(Ndim)/2
    
    M_log = Wad.kummer_log(a, c,k)
    Mplus_log = Wad.kummer_log(a + 1, c +1,k)
    dM_log = np.log((a/c)) + Mplus_log

    g = np.exp(dM_log - M_log)
#    kummer = 
#    print Ndim, k, r
    return g - r

def Newton_kappa(kappa0,Ndim,r, Ninter = 10):
    kappa = kappa0
    a = 0.5
    c = float(Ndim)/2
    for i in range(Ninter):
        
        M = np.exp(Wad.kummer_log(a, c,kappa))
        Mplus = np.exp(Wad.kummer_log(a + 1, c +1,kappa))
        dM = (a/c)*Mplus 
#        dM = (a - c)*Mbplus/c + M
        g = dM/M
#        print g
        dg =  (1 - c/kappa)*g + (a/kappa) - g*g
        
        kappa = kappa - (g - r)/dg
        
#        print kappa
    return kappa
    
def Newton_kappa_log(kappa0,Ndim,r, Ninter = 10):
    kappa = kappa0
    a = 0.5
    c = float(Ndim)/2
    for i in range(Ninter):
        
        M_log = Wad.kummer_log(a, c,kappa)
        Mplus_log = Wad.kummer_log(a + 1, c +1,kappa)
        dM_log = np.log((a/c)) + Mplus_log
    
        g = np.exp(dM_log - M_log)
#        print g
        dg =  (1 - c/kappa)*g + (a/kappa) - g*g
        
        kappa = kappa - (g - r)/dg
        
#        print kappa
    return kappa
    
def get_Watson_muKappa_ML(X, rk = None):
    # This function obtains both efficiently and checking the sign and that
    # If No rk specified it is just one 
    n,d = X.shape
    if(type(rk) == type(None)):
        rk = np.ones((n,1))
        
    Sk, D,V = get_eigenDV_ML(X, rk = rk)
    
    # Solve this thing
    d_max = np.argmax(D)
    d_min = np.argmin(D)
    mu_pos = V[:,d_max]
    mu_neg = V[:,d_min]
    ## We solve the positive and the negative situations and output the one with
    ## the highest likelihood ? 
    
    ### TODO: Warning it could happen that the eigenvalue with the lowest variance,
    ## The variance is so low that it is 0, so maybe we should not consider the negative
    # case when this happens ? What should be the tolerance ?
    
    eigenValue_pos = D[d_max]
    eigenValue_min = D[d_min]
    
    # If the negative mu variance is too small, we dont compute it.
    # TODO: maybe another safe measure for r = 1 ? in the positive case ?
    r_neg = np.dot(mu_neg.T,Sk).dot(mu_neg)
    r_pos = np.dot(mu_pos.T,Sk).dot(mu_pos)
    
#    print (r_neg, r_pos)
    
    # TODO: if r_neg -> 0 and r_pos -> 1 we are dealing with a fully degnerated cluster 
    # that is focused in one sample.
    
    tolerance =  1e-3
    if (r_neg < tolerance and r_pos > 1-tolerance):
        # Case where we have a degenerated cluster
#        print "Degenerated cluster"
        raise RuntimeError('Degenerated cluster focus in one sample. Percentage_samples = %f', "Degenerated_Cluster_Error",np.sum(rk)/n)
        
    elif (r_neg < tolerance and r_pos < 1-tolerance):
        # Case where the negative kappa case is very unilikely
        kappa_pos = get_Watson_kappa_ML(X, mu_pos,  Sk = Sk, rk = rk)
        kappa = kappa_pos
        mu = mu_pos
#        print "TODO: Warninig we do not try negative Kappa coz explained variance is too low !"
    elif (r_neg > tolerance and r_pos > 1-tolerance):
        # Case where the positive kappa case is very unilikely
    
        kappa_neg = get_Watson_kappa_ML(X, mu_neg,  Sk = Sk, rk = rk)
        kappa = kappa_neg
        mu = mu_neg
    else:
        # Case where both are possible.
#        print "We can consider negative kappa coz, the eingen value explains some variance"
    #    print d_max, d_min
    #    print D[d_max], D[d_min]
    #    print mu_pos
    #    print mu_neg
        
        kappa_pos = get_Watson_kappa_ML(X, mu_pos,  Sk = Sk, rk = rk)
        kappa_neg = get_Watson_kappa_ML(X, mu_neg,  Sk = Sk, rk = rk)
    
     
    #    likelihood_pos = np.sum(Watson_pdf_log(X.T,mu_pos,kappa_pos))
    #    likelihood_neg = np.sum(Watson_pdf_log(X.T,mu_neg,kappa_neg))
    
    #    likelihood_pos = np.sum(np.exp(Watson_pdf_log(X.T,mu_pos,kappa_pos))*rk.T)
    #    likelihood_neg = np.sum(np.exp(Watson_pdf_log(X.T,mu_neg,kappa_neg))*rk.T)
    
        # The maximum weighted likelihood estimator
        likelihood_pos = np.sum(Wad.Watson_pdf_log(X.T,mu_pos,kappa_pos)*rk.T)
        likelihood_neg = np.sum(Wad.Watson_pdf_log(X.T,mu_neg,kappa_neg)*rk.T)

     
    #    print likelihood_pos, likelihood_neg
        if (likelihood_pos > likelihood_neg):
            if (D[0] == D[1]):
                print "Warning: Eigenvalue1 = EigenValue2 in MLmean estimation"
            kappa = kappa_pos
            mu = mu_pos
        else:
            if (D[0] == D[1]):
                print "Warning: Eigenvalue1 = EigenValue2 in MLmean estimation"
            kappa = kappa_neg
            mu = mu_neg
    return mu, kappa

# This function obtains the correlation matrix and its eigenvector and values
def get_eigenDV_ML(X, rk = None):
    n,d = X.shape
    if(type(rk) == type(None)):
        rk = np.ones((n,1))
        
#    print (X*rk).shape
    Sk = np.dot(X.T,X*rk)   # Correlation
    Sk = Sk/np.sum(rk)            # Not really necesarry

    # Get eigenvalues to obtain the mu
    try:
        # Maybe rk is very small and this fails
        D,V = np.linalg.eig(Sk) # Obtain eigenvalues D and vectors V
    
        # Sometimes this gives also imaginaty partes 
        # TODO: Check what causes the imaginarity
        D = D.astype(float)
        V = V.astype(float)
    except np.linalg.linalg.LinAlgError:
        print "Sk has failed to have good parameters "
        print Sk
    
#    print D
    
    return Sk, D,V #  mu_pos, mu_neg

# This function obtains the positive and negative mus
def get_Watson_mus_ML(X, rk = None):
    Sk, D,V = get_eigenDV_ML(X, rk = rk)
    # Solve this thing
    d_max = np.argmax(D)
    d_min = np.argmin(D)
    mu_pos = V[:,d_max]
    mu_neg = V[:,d_min]
    return  mu_pos, mu_neg
    
def get_Watson_kappa_ML(X, mu,  Sk = None, rk = None):
    n,d = X.shape
    a = 0.5
    c = float(d)/2
    
    if (type(Sk) == type(None)):
        Sk = np.dot(X.T,rk*X)   # Correlation # We weight the samples by the cluster responsabilities r[:,k]
        Sk = Sk/(np.sum(rk))
    
    r = np.dot(mu.T,Sk).dot(mu)
    
    if (r == 1):
        r = r - 1e-10
        print "r = 1, percentage of samples: %f" %(100*np.sum(rk)/n)
#        raise RuntimeError("r = 1, we cannor continue. ercentage of samples %f" %(100*np.sum(rk)/n))
    # General aproximation
    BGG = (c*r -a)/(r*(1-r)) + r/(2*c*(1-r))
    # When r -> 1 
#    BGG = (c - a)/(1-r) + 1 - a + (a - 1)*(a-c-1)*(1-r)/(c-a)
    # TODO: In some examples this does not converge in the scipy method..

    kappa_max_warn = 1000
    if (np.abs(BGG) > kappa_max_warn):
#        warnings.warn("Initial guess of kappa = %f is too big.r: %f. Percentage of samples %f" %(BGG,r, 100*np.sum(rk)/n),UserWarning, stacklevel=2)
        pass
#        print "mu", np.sum(mu*mu), np.sum(mu), mu
 
#    try:
#    BGG_opt = newton(get_kappaNewton, BGG, args=([d,r],))
    BGG_opt = Newton_kappa_log(BGG,d,r,Ninter = 30)
#    BGG_opt = BGG
    if (np.abs(BGG_opt) > kappa_max_warn):
#            warnings.warn("Kappa after optimization %f is too big. Percentage of samples %f" %(BGG_opt, 100*np.sum(rk)/n),UserWarning, stacklevel=2)
        pass
    pass
    
#    except RuntimeError:
#        print "Oops! Newton did not converge"


    return BGG_opt
    

