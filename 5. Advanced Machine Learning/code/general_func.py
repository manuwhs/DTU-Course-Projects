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
import copy

def normalize_module(Xdata):
#    tol = 0.0000001
    # Expects a matrix (Nsamples, Ndim) and normalizes the values
#    print Xdata.shape
    Nsamples, Ndim = Xdata.shape
#    mean_of_time_instances = np.mean(Xdata, axis = 0).reshape(1,Ndim)
#    mean_of_channels = np.mean(Xdata, axis = 1).reshape(Nsamples,1)
#     Substract mean
#    Xdata = Xdata - 
    Xdata = Xdata  # - mean_of_channels# - mean_of_time_instances # - mean_of_time_instances# - mean_of_time_instances
    
#    print Xdata.shape
    # Normalice module
    Module = np.sqrt(np.sum(np.power(Xdata,2), axis = 1))
    Module = Module.reshape(Nsamples,1)
    # Check that the modulus is not 0
#    Xdata = Xdata[np.where(Module > tol)[0],:]
    Xdata = np.divide(Xdata,Module)
    
    return Xdata

def accuracy (Y,T):
    N_samples = len(Y)
    score = 0
    for i in range (N_samples):
#            print predicted[i], Y[i]
        if (Y[i] == T[i]):
            score += 1;
    return 100*score/float(N_samples)
    

def draw_HMM_indexes(pi, A, Nchains = 10, Nsamples = 30):
    # If Nsamples is a number then all the chains have the same length
    # If it is a list, then each one can have different length
    K = pi.size  # Number of clusters
    Chains_list = []
    
    Cluster_index = range(K)
    
    Nsamples = ul.fnp(Nsamples)
    if(Nsamples.size == 1):  # If we only have one sample 
        Nsamples = [int(Nsamples)]*Nchains
        
    for nc in range(Nchains):
        Chains_list.append([])
        sample_indx = np.random.choice(Cluster_index, 1, p = pi)
        Chains_list[nc].append(int(sample_indx))
        
        for isam in range(1,Nsamples[nc]):
            # Draw a sample according to the previous state
            sample_indx = np.random.choice(Cluster_index, 1, 
                                           p = A[sample_indx,:].flatten())
            Chains_list[nc].append(int(sample_indx))
    
    return Chains_list

def draw_HMM_samples(Chains_list, Samples_clusters):
    # We take the indexes of the chain and then draw samples from a pregenerated set
    # Samples_clusters is a list where every element is an array of sampled of the
    # i-th cluster
    
    K = len(Samples_clusters)

    Nchains = len(Chains_list)
    HMM_chains = [];
    
    counter_Clusters = np.zeros((K,1))
    for nc in range(Nchains):
        Nsamples = len(Chains_list[nc])
        HMM_chains.append([])
        
        for isam in range(0,Nsamples):
            K_index = Chains_list[nc][isam]
            Sam_index = int(counter_Clusters[K_index])
  
            sample = Samples_clusters[K_index][Sam_index,:]
            counter_Clusters[K_index] = counter_Clusters[K_index] +1
            HMM_chains[nc].append(sample)
    
        HMM_chains[nc] = np.array(HMM_chains[nc])
    return HMM_chains

def get_EM_data_from_HMM(HMM_list, Nchains_load = -1):
    # If Nchains_load = -1, it loads them all
    # Load first dataset
    k = 0 # For the initial
    for Xdata_chain in HMM_list:
        
        Xdata_chain = np.array(Xdata_chain)
#        gl.scatter_3D(Xdata_chain[:,0], Xdata_chain[:,1],Xdata_chain[:,2], nf = 0, na = 0)
    
    #    print Xdata_k.shape
        if (k == 0):
            Xdata = copy.deepcopy(Xdata_chain)
        else:
            Xdata = np.concatenate((Xdata, copy.deepcopy(Xdata_chain)), axis = 0)
        k += 1
        
        if k == Nchains_load: # Control how many samples we use
            break
        
    return Xdata
    