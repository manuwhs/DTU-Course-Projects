import Watson_distribution as Wad
import Watson_sampling as Was
import Watson_estimators as Wae
import general_func as gf

import numpy as np
import HMM_libfunc2 as HMMlf
#import HMM_libfunc2_old as HMMlf
import HMM_libfunc2_old as HMMlf_old

import copy

def HMM(data, I = 3,delta = 0.01,R = 30,
        pi_init = None, A_init = None, B_init = None,
        verbose = 0):
    # Input
    # S = Number of States
    # data = Data
    # alpha = minimu step for convergence (If negative it does not check it)
    # R = Max iterations
    # 
    # Output
    # pi = pi parameters
    # A =  A parameters
    # B =  B parameters
    # logL = complete log-likelihood of each iteration

    #########################################################
    # data = list of realizations, every realization has N samples !! 
    # It can vary from realization to realization

    N = len(data)       # Number of Realizations of the HMM
    D = data[0].shape[1]; # Dimension of multidimensial bernoulli
    T = data[0].shape[0]; # Number of samples of the HMM

    #********************************************************
    #*********** INITIALIZATION *****************************
    #********************************************************


    pi, A, B = HMMlf.init_HMM_params(D,I,pi_init, A_init, B_init, Kappa_max = 20)
    
    #Initialize log-likelihood to 0
    
    #*********************************************************
    #*********** ITERATIONS OF THE HMM ************************
    #*********************************************************
    logl = []   # List where we store the likelihoos
    B_list = []
    pi_list = []
    A_list = []
    
    pi_list.append(copy.deepcopy(pi))
    B_list.append(copy.deepcopy(B))
    A_list.append(copy.deepcopy(A))
    
    
    ## Get the samples in a line. We will need this later
    ######## THIS IS A CONCATENATED VERSION OF THE SAMPLES THAT WE NEED
    ## FOR SOME PART OF THE ALGO. But we also use data[0]. We take into account different chains.
    X = data[0]
    for n in range(1,N):
        X = np.concatenate((X, data[n]), axis = 0)
                
    for r in range(R):         # For every iteration of the EM algorithm
        if (verbose > 0):
            print "Iteration %i"%(r)
        
        ## We avoid 0s in A or pi...
        pi = pi + 1e-200
        A = A + 1e-200

        ##***********************************************
        ########## Stability Check #############
        ##***********************************************
        # Some models can kind of die if the clusters converge badly.
        # For example if we have too many clusters, one of them could go to an outlier,
        # And the resposibility of that point would be 100% from that cluster
        # And the parameters of the cluster cannot be properly computed.
        # It is up to the stability function to modify the parameters.
        # Maybe removing the cluster or reininitializing the cluster randomly
        
        B_new, pi_new, A_new, clusters_change = HMMlf.manage_clusters(X, B, pi, A)
        # If we reduced the 
        if (clusters_change):
            alpha = HMMlf.get_alfa_matrix_log(A,B,pi,data);

        #******************************************************  
        #*********** E Step ***********************************
        #******************************************************
        
#        print "pi paramters"
#        print pi
#        print "mus"
#        print B[0].T
#        print "kappas"
#        print B[1]
#        print "A"
#        print A
        
        # In this step we calculate the alfas, betas, gammas and fis matrices
        
        if (r == 0):
            ## ALPHA is recomputed ar the end to 
            alpha = HMMlf.get_alfa_matrix_log(A,B,pi,data);
            # Compute the initial incomplete-loglikelihood
            ll = HMMlf.get_HMM_Incomloglike(A,B,pi,data,alpha)
            logl.append(ll)
            if (verbose > 0):
                print "Initial loglikelihood: %f" % ll

        # Probabilities get vanishingly small as T -> 0
        # Maybe take logarithms ?
        
#        print alpha[0,0,:]  # I * N * T
#        print alpha[0,1,:]
#        print alpha[0].shape
            
        beta = HMMlf.get_beta_matrix_log(A,B,pi,data);
#        beta = HMMlf_old.get_beta_matrix_log( A,B,data);
#        print "beta diff %f" % np.sum(np.sum(np.abs(beta_2[0] - beta[0])))
        
        # Probabilities get vanishingly small as t -> 0
        # Maybe take logarithms ?
#        print beta.shape
#        print alpha[0,0,:]  # I * N * T
#        print beta[0][0,1]
#        print beta[0].shape
        
        gamma = HMMlf.get_gamma_matrix_log(alpha,beta );
#        gamma_2 = HMMlf_old.get_gamma_matrix_log(alpha,beta );
#        print "gamma diff %f" % np.sum(np.sum(np.abs(gamma_2[0] - gamma[0])))
        
#        print gamma[0].shape
#        print gamma [0,0,:]
        fi = HMMlf.get_fi_matrix_log(A,B, alpha,beta,data );
#        fi_2 = HMMlf_old.get_fi_matrix_log(A,B, alpha,beta,data );
#        print "fi diff %f" % np.sum(np.sum(np.abs(fi_2[0] - fi[0])))
        
#        print fi [0,0,:]
        
        #*****************************************************   
        #*********** M Step ***********************************
        #*****************************************************
        # In this step we calculate the next parameters of the HMM
        
#        print gamma[0]
        for n in range(N):
            gamma[n] = np.exp(gamma[n])
            fi[n] = np.exp(fi[n])
#        print "-----------------------------$"
#        print gamma[0]
        # Calculate new initial probabilities
        pi = HMMlf.get_pi(gamma)
#        print "pi"
#        print pi
        
        # Calculate transition probabilities A

        A = HMMlf.get_A(fi)
#        print "A"
#        print A
        # Calculate the paramters B
    
        B = HMMlf.get_B(X, gamma, B)
        
        ## We avoid 0s in A or pi...
        pi = pi + 1e-200
        A = A + 1e-200
        #********************************************************* 
        #****** Calculate Incomplete log-likelihood  *************
        #*********************************************************
    
        # Remember that the Incomplete log-likelihood could decrease with
        # the number of iterations at some stage since the EM algorith 
        # maximizes the Complete log-likelihood (They are different)
    
        # Calculate Incomplete log-likelihood with the Forward Algorithm
        alpha = HMMlf.get_alfa_matrix_log(A,B,pi,data);
#        alpha_2 = HMMlf_old.get_alfa_matrix_log(A,B,pi,data);
#        print "alpha diff %f" % np.sum(np.sum(np.abs(alpha_2[0] - alpha[0])))
        
        new_ll = HMMlf.get_HMM_Incomloglike(A,B,pi,data,alpha)
#        new_ll2 = HMMlf.get_HMM_Incomloglike_beta(A,B,pi,data)
#        print new_ll2 - new_ll
        
        if (verbose > 0):
            print "Loglkelihood: %f " % new_ll
        #***************************************************** 
        #****** Convergence Checking *************************
        #*****************************************************
 
        logl.append(new_ll)
        pi_list.append(copy.deepcopy(pi))
        B_list.append(copy.deepcopy(B))
        A_list.append(copy.deepcopy(A))
    
        if (clusters_change == 0): # If we did not have to delete clusters
            if(new_ll-ll <= delta): # Maybe abs
                break;
            else:
                ll = new_ll;
        else:
            ll = new_ll
#        print'  R:'
#        print r;

    return logl,B_list,pi_list, A_list


def run_several_HMM(data, I = 3,delta = 0.01,R = 30,
        pi_init = None, A_init = None, B_init = None, Ninit = 5,
        verbose = 0):
    print "HMM number 1/%i" % (Ninit)
# We make a first run of the HMM
    [logl,B_list,pi_list, A_list] = HMM(data,I,delta,R, pi_init, A_init, B_init, verbose = verbose);
    best_logl = logl;   
    best_pi = pi_list;     
    best_B = B_list;   
    best_A = A_list
    best_final_logll = logl[-1]
    
    if (Ninit > 1):
        for i in range(1,Ninit):
            print "HMM number %i/%i" % (i+1,Ninit)
            [logl,B_list,pi_list, A_list]  = HMM(data,I,delta,R, pi_init, A_init, B_init, verbose = verbose);
            if (logl[-1] > best_final_logll):
                best_logl = logl;   
                best_pi = pi_list;     
                best_B = B_list;   
                best_A = A_list
                best_final_logll = logl[-1]
    
    return [best_logl, best_B, best_pi, best_A]
