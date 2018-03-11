
import Watson_distribution as Wad
import Watson_sampling as Was
import Watson_estimators as Wae
import general_func as gf
import copy
import numpy as np

def sum_logs(log_vector, byRow = False):
    # This functions sums a vector of logarithms
    # alfa[i,n,t] = np.logaddexp(aux, alfa[i,n,t])

    if (byRow == False):
        # We just add all the components
        log_vector = np.array(log_vector).flatten()
        log_vector = np.sort(log_vector) # Sorted min to max
                 
        a0 = float(log_vector[-1])
        others = np.array(log_vector[:-1]).flatten()
        N = 1
    else:
        # log_vector = (Nrows, NvaluestoAdd)
        N = log_vector.shape[0]
        log_vector = np.sort(log_vector, axis = 1) # Sorted min to max
        a0 = log_vector[:,[-1]]
        others = log_vector[:,:-1]
        
        if (0):
            log_vector = np.array(log_vector)
            N = log_vector.shape[0]
            result = []
            for i in range(log_vector.shape[0]):
                result.append(sum_logs(log_vector[i,:]))
            
            result = np.array(result).flatten().reshape(N,1)
            
            return result

#    print np.exp(others - a0).shape
    if (byRow == True):
        caca = np.sum(np.exp(others - a0),axis = 1)
        caca = caca.reshape(N,1)
    else:
        try:
            caca = np.sum(np.exp(others - a0))
        except AttributeError:
            print a0
            print others
            print type(others)
            print log_vector
#    print caca.shape
    result = a0 + np.log(1 + caca)
    
    if (byRow == True):
#        print result.shape
        result = result.flatten().reshape(N,1)
    return result
    
def get_alfa_matrix_log( A,B,pi,data ):
    I = A.shape[0]
    N = len(data)
    D = data[0].shape[1]
    T = [len(x) for x in data]
    
    kappas = B[1]
    cp_logs = []
    
    for i in range(I):
        cp_logs.append(Wad.get_cp_log(D,kappas[:,i]))
    alfa = [];
    
    # Calculate first sample
    for n in range(N): # For every chain
        alfa.append(np.zeros((I,T[n])));
        for i in range(I):  # For every state
            alfa[n][i,0] = np.log( pi[:,i]) + Wad.Watson_pdf_log(data[n][0,:], B[0][:,i], B[1][:,i], cp_log = cp_logs[i]);

    # Calculate the rest of the alfas recursively
    for n in range(N):          # For every chain
        for t in range(1, T[n]):           # For every time instant
            aux_vec = np.log(A[:,:]) + alfa[n][:,[t-1]]
            alfa[n][:,[t]] = sum_logs(aux_vec.T,  byRow = True)
#            print sum_logs(aux_vec.T,  byRow = True).shape
#            print alfa[n][:,[t]].shape
            alfa[n][:,[t]] +=  Wad.Watson_K_pdf_log(data[n][[t],:].T, B[0][:,:], B[1][:,:], cp_logs).T
#                print np.log(Wad.Watson_pdf(data[n][t,:], B[0][:,i], B[1][:,i]))
#                print    np.log(Wad.Watson_pdf(data[n][t,:], B[0][:,i], B[1][:,i]))# alfa[i,n,t] 

#            for i in range(I):      # For every state
#                aux_vec = np.log(A[:,[i]]) + alfa[n][:,[t-1]]
#                alfa[n][i,t] = sum_logs(aux_vec)
#                alfa[n][i,t] =  Wad.Watson_pdf_log(data[n][[t],:].T, B[0][:,i], B[1][:,i], cp_log = cp_logs[i]) + alfa[n][i,t] ;
                
    return alfa
    
def  get_beta_matrix_log( A,B,pi,data ):
    I = A.shape[0]
    N = len(data)
    D = data[0].shape[1]
    T = [len(x) for x in data]
    
    kappas = B[1]
    cp_logs = []
    
    for i in range(I):
        cp_logs.append(Wad.get_cp_log(D,kappas[:,i]))
    
    beta = [];
    
    # Calculate the last sample
    for n in range(N): # For every chain
        beta.append(np.zeros((I,T[n])));
        Nsam, Nd = data[n].shape
#        pi_end = get_final_probabilities(pi,A,Nsam)
        
        for i in range(I):
            beta[n][i,-1] = 0 # np.log( pi_end[0,i])
#        print beta[n][:,-1]
#            beta[n][i,-1] = np.log( pi_end[0,i]) + Wad.Watson_pdf_log(data[n][-1,:], B[0][:,i], B[1][:,i], cp_log = cp_logs[i]);

#            aux_vec = []
#            for j in range(J):
#                aux_vec.append(pi_end[:,i])
#                np.log( pi_end[:,i]) + Wad.Watson_pdf_log(data[n][0,:], B[0][:,i], B[1][:,i], cp_log = cp_logs[i]);

    # Calculate the rest of the betas recursively
    for n in range(N):     # For every chain
        for t in range(T[n]-2,-1,-1):  # For every time instant backwards
            aux_vec = np.log(A[:,:]) +  beta[n][:,[t+1]].T + \
            Wad.Watson_K_pdf_log(data[n][[t+1],:].T, B[0][:,:], B[1][:,:],cp_logs)
                
            beta[n][:,[t]] = sum_logs(aux_vec, byRow = True)
                
    return beta
    
def  get_fi_matrix_log( A,B, alpha= None,beta= None,data = None):
    I = A.shape[0]
    N = len(data)
    D = data[0].shape[1]
    T = [len(x) for x in data]
    
    kappas = B[1]
    cp_logs = []
    
    for i in range(I):
        cp_logs.append(Wad.get_cp_log(D,kappas[:,i]))

    fi = []
    
    for n in range(N):
        fi.append(np.zeros((I,I,T[n]-1)))

        zurullo = np.log(A[:,:])
        zurullo = zurullo.reshape(zurullo.shape[0],zurullo.shape[1],1)
        zurullo = np.repeat(zurullo,T[n]-1,axis = 2)
        
        mierda1 = beta[n][:,1:] 
        mierda1 = mierda1.reshape(1,mierda1.shape[0],mierda1.shape[1])
        mierda1 = np.repeat(mierda1,I,axis = 0)
        
        mierda2 = alpha[n][:,:-1] 
        mierda2 = mierda2.reshape(mierda2.shape[0],1,mierda2.shape[1])
        mierda2 = np.repeat(mierda2,I,axis = 1)
        
        caca = Wad.Watson_K_pdf_log(data[n][1:,:].T, B[0][:,:], B[1][:,:], cp_logs).T
        caca = caca.reshape(1,caca.shape[0],caca.shape[1])
        caca = np.repeat(caca,I,axis = 0)
        
        fi[n][:,:,:] = zurullo + caca + mierda1 + mierda2
        
    for n in range(N):
        for t in range (0, T[n]-1):
            # Normalize to get the actual fi
            fi[n][:,:,t] = fi[n][:,:,t] - sum_logs(fi[n][:,:,t]);  

    return fi
    
def get_gamma_matrix_log( alpha,beta ):
    I = alpha[0].shape[0]
    N = len(alpha)
    T = [x.shape[1] for x in alpha]

    gamma = []
    
    for n in range(N):
        gamma.append(np.zeros((I,T[n])))
        for t in range (0, T[n]):
            gamma[n][:,t] = alpha[n][:,t] + beta[n][:,t];
    
    for n in range(N):
        for t in range(T[n]):
            #Normalize to get the actual gamma
            gamma[n][:,t] = gamma[n][:,t] - sum_logs(gamma[n][:,t]);  

    return gamma

def get_pi(gamma):
    
    # Calculate new initial probabilities
    N = len(gamma)
    I = gamma[0].shape[0]

    pi = np.zeros((1,I))
    N_gamma = []
    for n in range(N):
        N_gamma.append (np.sum(gamma[n][:,0]));
    N_gamma = np.sum(N_gamma)
    
    for i in range(I):
        aux = []
        
        for n in range(N):
            aux.append(gamma[n][i,0])

        N_i_gamma = np.sum(aux)
        pi[0,i] = N_i_gamma/N_gamma;
        
    pi = pi.reshape(1,pi.size)
    return pi

def get_A(fi):

# Calculate transition probabilities A

    I = fi[0].shape[0]
    N = len(fi)
    A = -np.ones((I,I))

    for i in range(I):
#        print range(I)
        E_i_fi = []
        # Calculate vector ai = [ai1 ai2 ... aiJ]  sum(ai) = 1
        for n in range(N): 
            E_i_fi.append(np.sum(np.sum(fi[n][i,:,:])))
        E_i_fi = np.sum(E_i_fi)
        
        for j in range(I):
            E_ij_fi = []
            for n in range(N): 
                E_ij_fi.append(np.sum(fi[n][i,j,:]))
            E_ij_fi = np.sum(E_ij_fi)
            A[i,j] = E_ij_fi/E_i_fi;
            
#        print "A"
#        print A
  
    return A
    
def init_HMM_params(D,I,pi_init = None, A_init = None, B_init = None, Kappa_max = 20):
    # Here we will initialize the  parameters of the HMM, that is,the initial
    # probabilities of the state "pi", the transition probabilities "A" and the
    # parameters of the probability functions "B"

    # Initial probabilities
    # We set the Initial probabilities with uniform discrete distribution, this
    # way, the a priori probability of any vector to belong to any component is
    # the same.
    if (type(pi_init) == type(None)): # If not given an initialization
        pi = np.ones((1,I));
        pi = pi*(1/float(I));
    else:
        pi = np.array(pi_init).reshape(1,I)
        
    # Transition probabilities "A"
    # We set the Transition probabilities with uniform discrete distribution, this
    # way, the a priori probability of going from a state i to a state j is
    # the same, no matter the j.

    if (type(A_init) == type(None)): # If not given an initialization
        A = np.ones((I,I));   #A(i,j) = aij = P(st = j | st-1 = i)  sum(A(i,:)) = 1
        for i in range(I):
            A[i,:] =  A[i,:]*(1/float(I));
    else:
        A = A_init
    # Parameters of the probability functions "B"
    # Give random values to the transit parameters. Since in this case, all
    # theta parameters theta(d,k), are the Expected value of a Bernoulli, we
    # asign values to them at random accoding to a uniform continuous
    # distribution in the support (0,1).
    
    if (type(B_init) == type(None)): # If not given an initialization

        mus = np.random.randn(D,I);
        mus = gf.normalize_module(mus.T).T
        kappas = np.random.uniform(-1,1,I) * Kappa_max
        kappas = kappas.reshape(1,I)

    else:
        mus = np.array(B_init[0]).reshape((D,I))
        kappas = np.array(B_init[1]).reshape((1,I))
    # We store the parameter of the clusters in B
    B = copy.deepcopy([mus, kappas])
    
    return pi, A, B
    
def get_B(X, gamma, B, Kappa_max = 1000):
    N = len(gamma)
    I = gamma[0].shape[0]
    D = X.shape[1]
    mus = B[0]
    
    kappas = np.zeros((1,I))
    new_mu = np.zeros((D,I))
    for i in range(I): # For every cluster
         # We compute the gamma normalization of the cluster
        # The contribution of every sample is weighted by gamma[i,n,t];
        # The total responsibility of the cluster for the samples is N_i_gamma

        rk = gamma[0][i,:]
        for n in range(1,N):
            rk = np.concatenate((rk, gamma[n][i,:]), axis = 0)
#            rk = np.array(rk)
        rk = rk.reshape(rk.size,1)
        
#        new_mu = Wae.get_Weighted_MLMean(rk,X)  
        
        try:
            new_mu, kappas[:,i] = Wae.get_Watson_muKappa_ML(X, rk)

    #        print new_mu
            
        except RuntimeError as err:
            error_type = err.args[1]
            print err.args[0] % err.args[2]
            print """We saturate kappa to %f. But in the next estimation the estimated kappa will also be as bad"""% (Kappa_max)
            
            # TODO: This could not work if the Kappa_max is still to high
            
            new_mu_pos, new_mu_neg = Wae.get_Watson_mus_ML(X, rk)
            if (B[1][0,i] >= 0):
                new_mu = new_mu_pos
            else:
                new_mu = new_mu_neg
            
            kappas[:,i] = Kappa_max # * kappas[:,k]/np.abs(kappas[:,k])
            
        signs = np.sum(np.sign(new_mu *  mus[:,i]))
#            print signs
        if (signs < 0):
            mus[:,i] = -new_mu
        else:
            mus[:,i] = new_mu
            
#        kappas[:,i] = Wae.get_Weighted_MLkappa(rk,mus[:,i],X)
    
    B = [mus, kappas]
    return B

def get_final_probabilities(pi,A,N):
    Af = A
    for i in range(N-1):
        Af = Af.dot(A)
#        print Af
    pif = pi.dot(Af)

    return pif
    
def get_HMM_Incomloglike(A,B,pi,data, alpha = []):

    N = len(data)
    I = pi.size
    # Check if we have been given alpha so we do not compute it
    if (len(alpha) == 0):
        alpha = get_alfa_matrix_log(A,B,pi,data)
    new_ll = 0

    for n in range(N):    # For every HMM sequence
        ## Generate probablilities of being at the state qt = j in the last point
        # of the chainfs  
        Nsam, Nd = data[n].shape
        pi_end = get_final_probabilities(pi,A,Nsam)
        all_val = []
#        print pi_end.shape
#        print np.log(pi_end)
#        print I
        for i in range(I):
            all_val.append(alpha[n][i,-1]) # + np.log(pi_end)[0,i]  + np.log(pi_end)[0,i]
#            print np.log(pi_end)[0,i]
#        print all_val
        new_ll = new_ll +  sum_logs(all_val)#  sum_logs(alpha[n][:,-1] + np.log(pi_end).T);
        
    return new_ll

def get_HMM_Incomloglike_beta(A,B,pi,data, beta = []):

    N = len(data)
    I = pi.size
    # Check if we have been given alpha so we do not compute it
    if (len(beta) == 0):
        beta = get_beta_matrix_log(A,B,pi,data)
    new_ll = 0
    cp_logs = []
    I = A.shape[0]
    N = len(data)
    D = data[0].shape[1]
    T = [len(x) for x in data]
    
    kappas = B[1]
    for i in range(I):
        cp_logs.append(Wad.get_cp_log(D,kappas[:,i]))

    for n in range(N):    # For every HMM sequence
        ## Generate probablilities of being at the state qt = j in the last point
        # of the chainfs  
        Nsam, Nd = data[n].shape
        all_val = []
#        print pi_end.shape
#        print np.log(pi_end)
#        print I

        for i in range(I):
            all_val.append(beta[n][i,0] + np.log(pi)[0,i] + Wad.Watson_pdf_log(data[n][0,:], B[0][:,i], B[1][:,i], cp_log = cp_logs[i]))
#            print np.log(pi_end)[0,i]
#        print all_val
        new_ll = new_ll +  sum_logs(all_val)#  sum_logs(alpha[n][:,-1] + np.log(pi_end).T);
        
    return new_ll
    
def get_errorRate(real, pred):
    Nfails = 0
    Ntotal = 0
    for i in range(len(real)):
        T = np.array(real[i]).size
        
        for t in range(T):
    #        print decoded[i][t]
    #        print  HMM_list[i][t]
            if (int(real[i][t]) != int( pred[i][t])):
                Nfails += 1
            Ntotal += 1
    Failure = 100 * float(Nfails)/Ntotal
    return Failure

def remove_A(A, k):
    # Remove one of the states from A, we delete from both dimensions
    # and renormalize
    A = np.delete(A, k, axis = 1)
    A = np.delete(A, k, axis = 0)
    Asum = np.sum(A, axis = 1)
    A = A/ Asum
    return A
    
def remove_cluster( B, pi, A, k):
    # This function removed the cluster k from the parameters
    kappa = B[1][0,k]
    B[0] = np.delete(B[0], k, axis = 1)
    B[1] = np.delete(B[1], k, axis = 1)
    pi = np.delete(pi, k, axis = 1)
    pi = pi / np.sum(pi)
    A = remove_A(A, k)
    
    print "$$$$$$$$$$$$ cluster %i removed" % (k)
    print " Its kappa was %f"%(kappa)
#    print A.shape
#    print pi.shape
#    print B[0].shape
#    print B[1].shape
    
    return B, pi, A
#    print theta[0].shape
def manage_clusters(X, B, pi, A, Kappa_max = 1000):
    kappas = B[1]
    K = kappas.shape[1]
    Nsam,D = X.shape
    clusters_change = 0
    # Truco del almendruco
    kummer_check = []
    for k in range(K):
        try:
            kummer_check.append(Wad.get_cp_log(D,kappas[:,k]))
        except RuntimeError as err:
            print "Error in Managing clusters"
            error_type = err.args[1]
            print err.args[0] % err.args[2]
            print """We saturate kappa to %f. But in the next estimation the estimated kappa will also be as bad"""% (Kappa_max)
            
            # TODO: This could not work if the Kappa_max is still to high
            kappas[:,k] = Kappa_max * kappas[:,k]/np.abs(kappas[:,k])
            clusters_change = 1
    
    # We go backwards in case we erase, not to affect indexes
    for k in range(K):
#        
#        if (kappas[0,K - 1 -k] > 1000):
#           B, pi, A = remove_cluster(B, pi, A,K - 1 -k)
        pass
#        if (kummer_check[K - 1 -k] == 0): # If we fucked up
#            theta, pi = remove_cluster(theta,pi,K - 1 -k)

    # Maybe condition on pi as well ? If pi is too small.

#    K = B[0].shape[1]
#    for k in range(K):
#        if (np.sum(A[:,K - 1 -k]) <  0.01):  # Less than 5 percent of the states go to this one
#            B, pi, A = remove_cluster(B, pi, A,K -1 -k)
#    
    return B, pi, A, clusters_change
    
def match_clusters(mus_1, kappas_1, mus_2, kappas_2):
    # The index of the new clusters do not have to match with the order
    # of our previous cluster so we will assign to each new cluster the index
    # that is most similar to the ones we had 
    
    # Initially we will chose the one with lower distance between centroids
    
    pass

def get_initial_HMM_params_from_EM(EM_params):
    pimix = EM_params[0]
    theta = EM_params[1]
    
    I = pimix.size
    
    pi_init = pimix
    B_init = theta
    A_init = np.repeat(pi_init, I, axis = 0)
    
    return pi_init, B_init, A_init
    