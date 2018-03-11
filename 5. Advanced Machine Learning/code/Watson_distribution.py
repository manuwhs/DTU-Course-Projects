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


def get_cp(Ndim, kappa):
    gammaValue = gamma(float(Ndim)/2)
    M = hyp1f1(0.5, float(Ndim)/2, kappa)   # Confluent hypergeometric function 1F1(a, b; x)
    cp = gammaValue / (np.power(2*np.pi, float(Ndim)/2)* M)
    return cp
    ## TODO: Make it one func
    #def get_K_cps_log(D,kappas):
    #    # This function will compute the constants for several clusters
    #    cp_logs = []
    #    K = kappas.shape[1]
    #    for k in range(K):
    #        cp_logs.append(difu.get_cp_log(D,kappas[:,k]))
def get_cp_log(Ndim, kappa):
    gammaValue_log = np.log(gamma(float(Ndim)/2))
    # Confluent hypergeometric function 1F1(a, b; x)
    M_log = kummer_log(0.5, float(Ndim)/2, kappa)   
    cp_log = gammaValue_log - (np.log(2*np.pi) *(float(Ndim)/2) + M_log)

    return cp_log
    
def check_Kummer(Ndim, kappa):
    # This functions checks if the Kummer function will go to inf
    # Returns 1 if Kummer is stable, 0 if unstable
    f = hyp1f1(0.5, float(Ndim)/2, kappa) 
    if (np.isinf(f) == False):
        return 1
    else:
        return 0
        
def kummer_log(a,b,x):
    ## First try using the funcion in the library.
    ## If it is 0 or inf then we try to use our own implementation with logs
    ## If it does not converge, then we return None !!

    f = hyp1f1(a,b,x)
    if (np.isinf(f) == True):
#        warnings.warn("hyp1f1() is 'inf', trying log version,  (a,b,x) = (%f,%f,%f)" %(a,b,x),UserWarning, stacklevel=2)
        f_log = kummer_own_log(a,b,x)
#        print f_log
        
    elif(f == 0):
#        warnings.warn("hyp1f1() is '0', trying log version, (a,b,x) = (%f,%f,%f)" %(a,b,x),UserWarning, stacklevel=2)
        raise RuntimeError('Kummer function is 0. Kappa = %f', "Kummer_is_0", x)
#        f_log = kummer_own_log(a,b,x)  # TODO: We cannot do negative x, the functions is in log
    else:
        f_log = np.log(f)
#        print (a,b,x)
#        print f_log
        
    f_log = float(f_log)
    return f_log

def kummer_own_log(a,b,x):
    # Default tolerance is tol = 1e-10.  Feel free to change this as needed.
    tol = 1e-10;
    log_tol = np.log(tol)
    # Estimates the value by summing powers of the generalized hypergeometric
    # series:
    #      sum(n=0-->Inf)[(a)_n*x^n/{(b)_n*n!}
    # until the specified tolerance is acheived.
    
    log_term = np.log(x) + np.log(a) - np.log(b)
#    print a,b,x
#    f_log =  HMMl.sum_logs([0, log_term])
    
    n = 1;
    an = a;
    bn = b;
    nmin = 5;
    
    terms_list = []
    
    terms_list.extend([0,log_term])
    d = 0
    while((n < nmin) or (log_term > log_tol)):
      # We increase the n in 10 by 10 reduce overheading of  while
      n = n + d;
#      print "puto n %i"%(n)
#      print f_log
      an = an + d;
      bn = bn + d;
      
      d = 1
#      term = (x*term*an)/(bn*n);
      log_term1 = np.log(x) + log_term  + np.log(an+d) - np.log(bn+d) - np.log(n+d)
      d += 1
      log_term2 = np.log(x) + log_term1  + np.log(an+d) - np.log(bn+d) - np.log(n+d)
      d += 1
      log_term3 = np.log(x) + log_term2  + np.log(an+d) - np.log(bn+d) - np.log(n+d)
      d += 1
      log_term4 = np.log(x) + log_term3  + np.log(an+d) - np.log(bn+d) - np.log(n+d)
      d += 1
      log_term = np.log(x) + log_term4  + np.log(an+d) - np.log(bn+d) - np.log(n+d)
  
      terms_list.extend([log_term1,log_term2,log_term3,log_term4,log_term] )
      
      if(n > 10000):  # We fucked up
#        print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4"
#        print " Not converged "
#        print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4"
        # If we could not compute it, we raise an error...
        raise RuntimeError('Kummer function not converged after 10000 iterations. Kappa = %f', "Kummer_is_inf",x)
    f_log = HMMl.sum_logs(terms_list);
#    print "f_log success %f " % f_log
#    print "-----------------------------------------"
#    print n
#    print "-----------------------------------------"
    return f_log

def Watson_pdf (alpha, mu, kappa, cp = None):
    # Watson pdf for a 
    # mu: [mu0 mu1 mu...] p-1 dimesion angles in radians
    # kappa: Dispersion value
    # alpha: Vector of angles that we want to know the probability
####
    # Just make sure the matrixes are aligned
    mu = np.array(mu)
    mu = mu.flatten().reshape(mu.size,1)
    
    alpha = np.array(alpha)
    alpha = alpha.reshape(mu.size,alpha.size/mu.size)
    
    Ndim = mu.size #+ 1  ??
    
    # If we indicate cp, we do not compute it
    if (cp == None):
        cp = get_cp(Ndim, kappa)
#        print "GRGRGR"
#    print np.dot(mu.T, alpha)
    
    aux1 = np.dot(mu.T, alpha)
#    aux2 = 0
#    for i in range(mu.size):
#        aux2 = aux2 + mu[i]*alpha[i]
#    print alpha
    
#    if (kappa < 0):
#        print "Warning: Kappa < 0"
        
    pdf = cp * np.exp(kappa * np.power(aux1,2))
    
    if (pdf.size == 1): # Turn it into a single number if appropiate
        pdf = float(pdf)
    return pdf
    
      
def Watson_pdf_log (alpha, mu, kappa, cp_log = None):
    # Compute this in case that the probability is too high or low for just one sample
    # cp is ok, we can calculate normaly and then put it log
    # cp goes to very high if low dimensions and high kappa
    
    # If we indicate cp_log  we do not compute it.
    
    # Watson pdf for a 
    # mu: [mu0 mu1 mu...] p-1 dimesion angles in radians
    # kappa: Dispersion value
    # alpha: Vector of angles that we want to know the probability
####
    # Just make sure the matrixes are aligned

    mu = np.array(mu)
    mu = mu.flatten().reshape(mu.size,1)
    
    alpha = np.array(alpha)
    alpha = alpha.reshape(mu.size,alpha.size/mu.size)
    
    Ndim = mu.size #+ 1  ??
    
    if (type(cp_log) == type(None)):
        cp_log = get_cp_log(Ndim, kappa)
        
#    print np.dot(mu.T, alpha)
    
    aux1 = np.dot(mu.T, alpha)
#    aux2 = 0
#    for i in range(mu.size):
#        aux2 = aux2 + mu[i]*alpha[i]
#    print alpha

    log_pdf = cp_log + (kappa * np.power(aux1,2))
    
    if (log_pdf.size == 1): # Turn it into a single number if appropiate
        log_pdf = float(log_pdf)
    return log_pdf


def Watson_K_pdf_log (alpha, mus, kappas, cps_log = None):
    # Extension of Watson_pdf_log in which we also accept several clusters
    # We have to be more restrict in this case and the parameters must be:
    # alpha(D,Nsamples)  mu(D,K) kappa(K) cp_log(K)
    # The result is (Nsamples, K)

    Ndim, Nsam = alpha.shape
    Ndim2, K = mus.shape
    
    if (type(cps_log) == type(None)):
        cps_log = get_cp_log(Ndim, kappas)
        
    kappas = np.array(kappas)
    kappas = kappas.reshape(kappas.size,1)
    cps_log = np.array(cps_log)
    cps_log = cps_log.reshape(cps_log.size,1)
    
    aux1 = np.dot(mus.T, alpha)
#    aux2 = 0
#    for i in range(mu.size):
#        aux2 = aux2 + mu[i]*alpha[i]
#    print alpha
    log_pdf = cps_log + (kappas * np.power(aux1,2))
    return log_pdf.T
    