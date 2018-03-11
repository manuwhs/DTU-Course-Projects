
from scipy.special import hyp1f1
from scipy.special import gamma
import numpy as np
import copy 
import HMM_libfunc2 as HMMlf
import EM_libfunc as EMlf
import time
    
import Watson_distribution as Wad
import Watson_sampling as Was
import Watson_estimators as Wae
import general_func as gf

def EM(X, K = 3, delta = 0.01, T = 30, 
       pi_init = None, theta_init = None, verbose = 0):
    #     Input
    #     K = Cluster Number
    #     data = Data
    #     alpha = minimu step for convergence
    #     T = Max iterations
    #     pi_init = Optional initilization of pi
    #     theta_init = Optional initilization of theta

    #     Output
    #     pimix = pi parameters
    #     theta =  theta parameters
    #     logL = complete log-likelihood of each iteration
    
    # function [pimix,theta,logl] = EM(K,data,alpha,T)
    
    if (type(X) == type([])): # If we are given a list we joing the elements
                              # We expect it to be filled with matrix (Nsam_i, Ndim_j)
        X = np.concatenate(X, axis = 0)
#        print X.shape
        
    N = X.shape[0] # Number of IDD samples
    D = X.shape[1] # Dimension of dimensions of the data
    
    # theta:  Matrix whose k-th column is the D-dimensional theta vector of
    # parameters for the k-th component theta(:,k) = [theta_k1,..., theta_kD]
    
    #********************************************************
    #*********** INITIALIZATION *****************************
    #********************************************************
    
    pimix, theta = EMlf.init_EM_params(D,K, pi_init, theta_init, Kappa_max = 50)
    #*********************************************************
    #*********** ITERATIONS OF THE EM ************************
    #*********************************************************
    logl = []   # List where we store the likelihoos
    theta_list = [] # Kappas
    pimix_list = []
    
    # Calculate initial log-likelihood
#    ll = EMlf.get_EM_Incomloglike_log(theta,pimix,X)
    ll = -1e500  # -inf
#    print "Initial loglk: %f"%(ll)
    
#    logl.append(ll)
    theta_list.append(copy.deepcopy(theta))
    pimix_list.append(copy.deepcopy(pimix))


    ##***********************************************
    ########## Stability Check #############
    ##***********************************************
    # Some models can kind of die if the clusters converge badly.
    # For example if we have too many clusters, one of them could go to an outlier,
    # And the resposibility of that point would be 100% from that cluster
    # And the parameters of the cluster cannot be properly computed.
    # It is up to the stability function to modify the parameters.
    # Maybe removing the cluster or reininitializing the cluster randomly
    
    theta_new, pimix_new, clusters_change = EMlf.manage_clusters(X, theta, pimix)
    # If we reduced the 


    for t in range(T):    #For every iteration of the EM algorithm
        # T is the maximum number of iterations, if we stop before due to convergence
        # we will break the loop
        if (verbose > 0):
            print "Iteration %i"%t
        
        ## We avoid 0s pimix...
        pimix = pimix + 1e-200

        # TODO Check if we have removed clusters so that we do not stop the next iteration
        # Because probably the system will have less likelihood
        #******************************************************  
        #*********** E Step ***********************************
        #******************************************************
        
        # In this step we calculate the responsibility of each sample i to each
        # component k. This gives a measure of how likeky is sample i, to
        # belong to the component K.
        
#        t0 = time.time()

        ## OPTIMIZATION, we join the r_log obtaining and the ll, they share most of code
#        r_log = EMlf.get_responsabilityMatrix_log(X,theta,pimix)
#        ll = EMlf.get_EM_Incomloglike_log(theta,pimix,X)
        
        r_log, new_ll = EMlf.get_r_and_ll(X,theta,pimix)
#        t1 = time.time()
#        print "Time Responsability: %f"%(t1-t0)
        
        r = np.exp(r_log)
        ## TODO: Check stability here given r
        
#        r = EMlf.get_responsabilityMatrix(X,theta,pimix)
        #*****************************************************   
        #*********** M Step ***********************************
        #*****************************************************
    
        # In this step we calculate the next parameters of the mixture mdoel
        #Calculate new pimix and update

        pimix = EMlf.get_pimix(r)

#        print "pimix"
#        print pimix
        
        # Calculate new thetas parameters and update

        theta = EMlf.get_theta(r, X , theta)
#        print 'mus:'
#        print mus
#        
#        print "kappas"
#        print kappas
        
        theta_new, pimix_new, clusters_change = EMlf.manage_clusters(X, theta, pimix)
#        print clusters_change
            
        #********************************************************* 
        #****** Calculate Incomplete log-likelihood  *************
        #*********************************************************
    
        # Remember that the Incomplete log-likelihood could decrease with
        # the number of iterations at some stage since the EM algorith 
        # maximizes the Complete log-likelihood (They are different)
    
#        Now this is merged into a single function in Responsability obtaining
#        We actually do one iteration more !! We detect it one iteration delayed
#        new_ll = EMlf.get_EM_Incomloglike_log(theta,pimix,X)
        if (verbose > 0):        
            print "Loglk: %f"%(new_ll)

#        print'  new_ll:'
#        print new_ll
        
#        logl[t] = new_ll;
    
        #***************************************************** 
        #****** Convergence Checking *************************
        #*****************************************************
        
        logl.append(new_ll) # This is for the previous one
        theta_list.append(copy.deepcopy(theta))
        pimix_list.append(copy.deepcopy(pimix))
    
        if (t == T-1):  # If we are done with this
            new_ll = EMlf.get_EM_Incomloglike_log(theta,pimix,X)
            logl.append(new_ll)
            break
        
        if (clusters_change == 0): # If we did not have to delete clusters
            if(np.abs(new_ll-ll) <= delta):

                # Compute the last Loglikelihood
                new_ll = EMlf.get_EM_Incomloglike_log(theta,pimix,X)
                logl.append(new_ll)
                break;
            else:
                ll = new_ll;
        else:
            ll = new_ll
#        print'  R:'
#        print r;
    
    print "Final ll: %f"% (logl[-1])
    return logl,theta_list,pimix_list

# Runs the ME N times and chooses the realization with the least loglihood

def run_several_EM(X, K = 3, delta = 0.01, T = 30, 
       pi_init = None, theta_init = None, Ninit = 5, verbose = 0):

    if (type(X) == type([])): # If we are given a list we joing the elements
                              # We expect it to be filled with matrix (Nsam_i, Ndim_j)
        X = np.concatenate(X, axis = 0)
#        print X.shape
    print "EM number 1/%i" % Ninit
# We make a first run of the HMM
    [logl,theta_list,pimix_list] = EM(X,K,delta,T, pi_init,theta_init, verbose = verbose);
    best_logl = logl;   
    best_pimix = pimix_list;     
    best_theta = theta_list;     
    best_final_logll = logl[-1]
    
    if (Ninit > 1):
        for i in range(1,Ninit):
            print "EM number %i/%i" % (i+1,Ninit)
            [logl,theta_list,pimix_list] = EM(X,K,delta,T, pi_init,theta_init, verbose = verbose)
            
            if (logl[-1] > best_final_logll):
                best_logl = logl;   
                best_pimix = pimix_list;     
                best_theta = theta_list;     
                best_final_logll = logl[-1]
    
    return [best_logl, best_theta, best_pimix]
