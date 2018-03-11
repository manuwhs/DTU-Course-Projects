
# Official libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Own libraries
import import_folders
from graph_lib import gl
import sampler_lib as sl
import EM_libfunc as EMlf
import EM_lib as EMl
import copy
import pickle_lib as pkl

import Watson_distribution as Wad
import Watson_sampling as Was
import Watson_estimators as Wae
import general_func as gf

plt.close("all")

mu_caca = np.ones((1000,1))

################################################################
######## Load and combine 3 sets ##############################
###############################################################

EM_data = 1
if (EM_data):
    K = 3
    #gl.scatter_3D([0,1,1,1,1,-1,-1,-1,-1], [0,1,1,-1,-1,1,1,-1,-1],[0,1,-1,1,-1,1,-1,1,-1], nf = 1, na = 0)
    gl.scatter_3D(0, 0,0, nf = 1, na = 0)
    kflag = 0
    for k in range(1,K+1):
#        folder = "./EM_data/"
        folder = "./test_data/"
        
        filedir = folder + "Wdata_"+ str(k)+".csv"
        Xdata_k = np.array(pd.read_csv(filedir, sep = ",", header = None))
        Xdata_k = Xdata_k[:1000,:]
        print Xdata_k.shape
    #    Xdata_param = pkl.load_pickle( folder + "Wdata_"+ str(k)+".pkl",1)
    #    mu = Xdata_param[0]
    #    kappa = Xdata_param[1]
        
    #    print "Real: ", mu,kappa
        # Generate and plot the data
        gl.scatter_3D(Xdata_k[:,0], Xdata_k[:,1],Xdata_k[:,2], nf = 0, na = 0)
        
        mu_est2, kappa_est2 = Wae.get_Watson_muKappa_ML(Xdata_k)
        print "ReEstimated: ", mu_est2,kappa_est2
        
        if (kflag == 0):
            Xdata = copy.deepcopy(Xdata_k)
            kflag = 1
        else:
            Xdata = np.concatenate((Xdata, copy.deepcopy(Xdata_k)), axis = 0)

################################################################
######## Or Load the same data as for HMM 3 sets ###############
###############################################################

HMM_data = 1
if (HMM_data):
    folder = "./HMM_data/"
    HMM_list = pkl.load_pickle(folder +"HMM_datapoints.pkl",1)
    #gl.scatter_3D([0,1,1,1,1,-1,-1,-1,-1], [0,1,1,-1,-1,1,1,-1,-1],[0,1,-1,1,-1,1,-1,1,-1], nf = 1, na = 0)
#    gl.scatter_3D(0, 0,0, nf = 1, na = 0)
    k = 0 # For the initial
    for Xdata_chain in HMM_list:
    
        Xdata_chain = np.array(Xdata_chain)
    #    print Xdata_k.shape
#        gl.scatter_3D(Xdata_chain[:,0], Xdata_chain[:,1],Xdata_chain[:,2], nf = 0, na = 0, color = "k")
        if (k == 0):
            Xdata = copy.deepcopy(Xdata_chain)
        else:
            Xdata = np.concatenate((Xdata, copy.deepcopy(Xdata_chain)), axis = 0)
        k = 1

################################################################
######## Perform the EM !! ###############
###############################################################
perform_EM = 1
if (perform_EM):
    K = 3
    D = Xdata.shape[1]
    
    pi_init = np.ones((1,K));
    pi_init = pi_init*(1/float(K));
    
    mus_init  = np.random.randn(D,K);
    mus_init  = gf.normalize_module(mus_init.T).T
    
    kappas_init = np.random.uniform(-1,1,K) * 10
    kappas_init = kappas_init.reshape(1,K)
    theta_init = [mus_init , kappas_init ]
    
    ## RUN only one !!
    logl,theta_list,pimix_list = EMl.EM(Xdata, K = K, delta = 0.1, T = 100,
                                    pi_init = pi_init, theta_init = theta_init)
    
    ## RUN several !!
    
    #logl,theta_list,pimix_list = EMl.run_several_EM(Xdata, K = K, delta = 0.1, T = 30,
    #                            Ninit = 5)
    
    
    mus_list = []
    kappas_list =[]
    for theta in theta_list:
        mus_list.append(theta[0])
        kappas_list.append(theta[1]) 
    
    print "kappas"
    print kappas_list[-1]
    print "mus"
    print mus_list[-1]
    print "pimix"
    print pimix_list[-1]
    
    plot_evolution = 1
    if (plot_evolution):
        # Only doable if the clusters dont die
        mus_array = np.array(mus_list) # (Nit,Ndim,Nk)
        Nit,Ndim,Nk = mus_array.shape
        for k in range(Nk):
            gl.scatter_3D(mus_array[:,0,k], mus_array[:,1,k],mus_array[:,2,k], nf = 0, na = 0, join_points = "yes", alpha = 0.1)
        
gl.plot(range(1,np.array(logl).flatten()[1:].size +1),np.array(logl).flatten()[1:], 
        legend = ["EM LogLikelihood"], 
labels = ["Convergence of LL with generated data","Iterations","LL"], 
lw = 4,
        fontsize = 25,   # The font for the labels in the title
        fontsizeL = 30,  # The font for the labels in the legeng
        fontsizeA = 20)

#gl.plot(range(1,np.array(caca).flatten()[1:].size +1),np.array(caca).flatten()[1:], nf = 0,
#        legend = ["HMM LogLikelihood"], 
#labels = ["Convergence of LL with generated data","Iterations","LL"], 
#lw = 4,
#        fontsize = 25,   # The font for the labels in the title
#        fontsizeL = 30,  # The font for the labels in the legeng
#        fontsizeA = 20)
#        
#Sel = EMlf.EMdecode(Xdata,[mus_list[-1], kappas_list[-1]], pimix_list[-1])

