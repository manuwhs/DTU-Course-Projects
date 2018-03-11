
# Official libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Own libraries
import import_folders
from graph_lib import gl

import sampler_lib as sl
import EM_lib as EMl
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
Chains_list = pkl.load_pickle(folder +"HMM_labels.pkl",1)
HMM_list = pkl.load_pickle(folder +"HMM_datapoints.pkl",1)
params = pkl.load_pickle(folder +"HMM_param.pkl",1)
pi = params[0]
A = params[1]
pi_end = HMMlf.get_final_probabilities(pi,A,20)
#print pi_end

print "Real pi"
print pi
print "Real A"
print A

gl.scatter_3D(0, 0, 0, nf = 1, na = 0)
for XdataChain in HMM_list:
    gl.scatter_3D(XdataChain[:,0], XdataChain[:,1],XdataChain[:,2], nf = 0, na = 0)


################################################################
######## Initialization of the parameters ##############################
###############################################################


I = 3
D = HMM_list[0].shape[1]

init_with_EM = 0

if (init_with_EM == 0):
    pi_init =  np.ones((1,I));
    pi_init = pi_init*(1/float(I));
    A_init = np.ones((I,I));   #A(i,j) = aij = P(st = j | st-1 = i)  sum(A(i,:)) = 1
    for i in range(I):
        A_init[i,:] =  A_init[i,:]*(1/float(I));
    
    mus_init = np.random.randn(D,I);
    mus_init = gf.normalize_module(mus_init.T).T
    kappas_init = np.random.uniform(-1,1,I) * 10
    B_init = [mus_init, kappas_init]

elif (init_with_EM):
    X = HMM_list[0]
    N = len(HMM_list)
    for n in range(1,N):
        X = np.concatenate((X, HMM_list[n]), axis = 0)
            
    logl,theta_list,pimix_list = EMl.EM(X, K = I, delta = 0.1, T = 30)
    
    pi_init = pimix_list[-1]
    B_init = theta_list[-1]
    A_init = np.repeat(pi_init, I, axis = 0)

## ############################## Run 1 ###############################
perform_HMM = 0
if (perform_HMM):
    logl,B_list,pi_list, A_list = \
        HMMl.HMM(data = HMM_list,I = I,delta = 0.01, R = 20
                 ,pi_init = pi_init, A_init = A_init, B_init = B_init)
    
    ## RUN several !!
    #logl,B_list,pi_list, A_list = \
    #    HMMl.run_several_HMM(data = HMM_list,I = I,delta = 0.01, R = 20,
    #             Ninit = 5)
    
    mus_list = []
    kappas_list =[]
    for B in B_list:
        mus_list.append(B[0])
        kappas_list.append(B[1]) 
    
    print "kappas"
    print kappas_list[-1]
    print "mus"
    print mus_list[-1]
    print "pi"
    print pi_list[-1]
    print "A"
    print A_list[-1]
    print "Likelihood"
    print logl[-1]
    
    plot_evolution = 1
    # Plot evolution of the clusters
    if (plot_evolution):
        # Only doable if the clusters dont die
        mus_array = np.array(mus_list) # (Nit,Ndim,Nk)
        Nit,Ndim,Nk = mus_array.shape
        for k in range(Nk):
            gl.scatter_3D(mus_array[:,0,k], mus_array[:,1,k],mus_array[:,2,k], nf = 0, na = 0, join_points = "yes")
    
gl.plot([],np.array(logl).flatten()[1:], 
        legend = ["HMM LogLikelihood"], 
labels = ["Convergence of LL with generated data","Iterations","LL"], 
lw = 4,
        fontsize = 25,   # The font for the labels in the title
        fontsizeL = 30,  # The font for the labels in the legeng
        fontsizeA = 20)

#caca = copy.deepcopy(logl)

decode_shit = 0
if (decode_shit == 1):
    decoded_SbS = decl.SbS_decoder(data = HMM_list,
                               A = A_list[-1], B = [mus_list[-1], kappas_list[-1]],
                               pi = pi_list[-1])
    
    decoded_ML = decl.MLViter_decoder(data = HMM_list,
                               A = A_list[-1], B = [mus_list[-1], kappas_list[-1]],
                               pi = pi_list[-1])
    
    decoded_MAP = decl.MAPViter_decoder(data = HMM_list,
                               A = A_list[-1], B = [mus_list[-1], kappas_list[-1]],
                               pi = pi_list[-1])
    
    Failure_SbS = HMMlf.get_errorRate(Chains_list, decoded_SbS)                          
    Failure_ML =  HMMlf.get_errorRate(Chains_list, decoded_ML)      
    Failure_MAP =  HMMlf.get_errorRate(Chains_list, decoded_MAP)      
    
    new_ll = HMMlf.get_HMM_Incomloglike(A = A_list[-1],
                                        B = [mus_list[-1], kappas_list[-1]],
                                        pi = pi_list[-1],data = HMM_list)

    print "likelihood = %f"%(new_ll)
    print "Failure SbS %f"%(Failure_SbS)
    print "Failure ML %f"%(Failure_ML)
    print "Failure MAP %f"%(Failure_MAP)
