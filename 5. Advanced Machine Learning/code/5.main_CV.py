
# Official libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import copy
# Own libraries
import import_folders
from graph_lib import gl

import sampler_lib as sl
import EM_lib as EMl
import EM_libfunc as EMlf
import HMM_lib as HMMl
import HMM_libfunc2 as HMMlf
import decoder_lib as decl
import pickle_lib as pkl
import general_func as gf
plt.close("all")

################################################################
######## Load and combine 3 sets ##############################
###############################################################

folder = "./HMM_data/"
HMM_list = pkl.load_pickle(folder +"HMM_datapoints.pkl",1)
HMM_list2 = pkl.load_pickle(folder +"HMM2_datapoints.pkl",1)

for i in range (len(HMM_list)):
    chain = HMM_list[i]
    Nsamples, Ndim = HMM_list[i].shape
    chain_noise = np.random.rand(Nsamples, Ndim)/10
    HMM_list[i] = gf.normalize_module(chain + chain_noise)

for i in range (len(HMM_list2)):
    chain = HMM_list2[i]
    Nsamples, Ndim = HMM_list2[i].shape
    chain_noise = np.random.rand(Nsamples, Ndim)/10
    HMM_list2[i] = gf.normalize_module(chain + chain_noise)
    
D = HMM_list[0].shape[1]

## Examples Crosvalidation fot HMM
CV_HMM = 0
if (CV_HMM):
    
    ## Run 1
    States = [1,2,3,4,5,6,7]
    final_logl_tr = 0
    final_logl_val = 0

    logl_tr = []
    logl_val = []
    # For X in trance
    for I in States:
        logl,B_list,pi_list, A_list = \
            HMMl.run_several_HMM(data = HMM_list,I = I,delta = 0.01, R = 20
                    ,Ninit = 10)
       
        A = A_list[-1]
        pi = pi_list[-1]
        B = B_list[-1]
        ## Compute the likelihoods for train and test
        new_ll = HMMlf.get_HMM_Incomloglike(A,B,pi ,data = HMM_list)
        logl_tr.append(copy.deepcopy(new_ll))
        new_ll = HMMlf.get_HMM_Incomloglike(A,B,pi ,data = HMM_list2)
        logl_val.append(copy.deepcopy(new_ll))
    logl_tr = np.array(logl_tr)
    logl_val = np.array(logl_val)
    final_logl_tr += logl_tr
    final_logl_val += logl_val

    logl_tr = []
    logl_val = []

    for I in States:
        logl,B_list,pi_list, A_list = \
            HMMl.run_several_HMM(data = HMM_list2,I = I,delta = 0.01, R = 20
                    ,Ninit = 5)
       
        A = A_list[-1]
        pi = pi_list[-1]
        B = B_list[-1]
        ## Compute the likelihoods for train and test
        new_ll = HMMlf.get_HMM_Incomloglike(A,B,pi ,data = HMM_list2)
        logl_tr.append(copy.deepcopy(new_ll))
        new_ll = HMMlf.get_HMM_Incomloglike(A,B,pi ,data = HMM_list)
        logl_val.append(copy.deepcopy(new_ll))
    logl_tr = np.array(logl_tr)
    logl_val = np.array(logl_val)
    final_logl_tr += logl_tr
    final_logl_val += logl_val

    final_logl_tr = final_logl_tr/2
    final_logl_val = final_logl_val/2
    
#    gl.plot(States, final_logl_tr, legend = ["tr"], labels = ["HMM","States","loglike"])
#    gl.plot(States, final_logl_val, nf = 0, legend = ["Val"])
    
    gl.plot(States,final_logl_tr, 
            legend = ["Train LL"], 
    labels = ["Validation of Number of clusters HMM with LL","Number of clusters (K)","LL"], 
    lw = 4,
            fontsize = 25,   # The font for the labels in the title
            fontsizeL = 30,  # The font for the labels in the legeng
            fontsizeA = 20)
    
    gl.plot(States,final_logl_val, nf = 0,
            legend = ["Validation LL"], 
    lw = 4,
            fontsize = 25,   # The font for the labels in the title
            fontsizeL = 30,  # The font for the labels in the legeng
            fontsizeA = 20)

#caca1 = copy.deepcopy(final_logl_tr)
#caca2 = copy.deepcopy(final_logl_val)
CV_EM = 1
if (CV_EM):
    ## Run 1
    Klusters = [1,2,3,4,5,6,7] # range(1,8) # 3,4,5,6,10,10,12,15
    final_logl_tr = 0
    final_logl_val = 0
    
    Xdata = gf.get_EM_data_from_HMM(HMM_list)
    Xdata2 = gf.get_EM_data_from_HMM(HMM_list2)

    logl_tr = []
    logl_val = []
    # For X in trance
    for K in Klusters:
        logl,theta_list,pimix_list = EMl.run_several_EM(Xdata, K = K, delta = 0.1, T = 50,
                                    Ninit = 10)     
        theta = theta_list[-1]
        pimix = pimix_list[-1]
        
        ## Compute the likelihoods for train and test
        new_ll = EMlf.get_EM_Incomloglike_log(theta,pimix = pimix_list[-1],X = Xdata)
        logl_tr.append(copy.deepcopy(new_ll))
        new_ll = EMlf.get_EM_Incomloglike_log(theta,pimix = pimix_list[-1],X = Xdata2)
        logl_val.append(copy.deepcopy(new_ll))
        
    logl_tr = np.array(logl_tr)
    logl_val = np.array(logl_val)
    final_logl_tr += logl_tr
    final_logl_val += logl_val

    logl_tr = []
    logl_val = []
    for K in Klusters:
        logl,theta_list,pimix_list = EMl.run_several_EM(Xdata2, K = K, delta = 0.1, T = 50,
                                    Ninit = 5)     
        theta = theta_list[-1]
        pimix = pimix_list[-1]
        
        ## Compute the likelihoods for train and test
        new_ll = EMlf.get_EM_Incomloglike_log(theta,pimix = pimix_list[-1],X = Xdata2)
        logl_tr.append(copy.deepcopy(new_ll))
        new_ll = EMlf.get_EM_Incomloglike_log(theta,pimix = pimix_list[-1],X = Xdata)
        logl_val.append(copy.deepcopy(new_ll))
        
    logl_tr = np.array(logl_tr)
    logl_val = np.array(logl_val)
    final_logl_tr += logl_tr
    final_logl_val += logl_val
    
    final_logl_tr = final_logl_tr/2
    final_logl_val = final_logl_val/2
    
    
    gl.plot(Klusters,final_logl_tr, 
            legend = ["Train EM"], 
    labels = ["Validation of Number of clusters with LL","Number of clusters (K)","LL"], 
    lw = 4,
            fontsize = 25,   # The font for the labels in the title
            fontsizeL = 30,  # The font for the labels in the legeng
            fontsizeA = 20)
    
    gl.plot(Klusters,final_logl_val, nf = 0,
            legend = ["Validation EM"], 
    lw = 4,
            fontsize = 25,   # The font for the labels in the title
            fontsizeL = 30,  # The font for the labels in the legeng
            fontsizeA = 20)


#    gl.plot(States,caca1, nf = 0,
#            legend = ["Train HMM"], 
#    lw = 4,
#            fontsize = 25,   # The font for the labels in the title
#            fontsizeL = 30,  # The font for the labels in the legeng
#            fontsizeA = 20)
#            
#    gl.plot(States,caca2, nf = 0,
#            legend = ["Validation HMM"], 
#    lw = 4,
#            fontsize = 25,   # The font for the labels in the title
#            fontsizeL = 30,  # The font for the labels in the legeng
#            fontsizeA = 20)
#            
#    gl.plot(Klusters, final_logl_tr, legend = ["tr"], labels = ["EM","K","loglike"])
#    gl.plot(Klusters, final_logl_val, nf = 0, legend = ["Val"])




# View final kappas theta_list[-1][1]