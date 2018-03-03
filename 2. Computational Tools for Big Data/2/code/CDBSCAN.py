import numpy as np
import copy 
import scipy.sparse as sps
import gc as gc

def remove_from_list(L,RL):
    # This function removes the elements RL from list L
    Lc = copy.copy(L)  # Create a copy of L
    for rl in RL:  # For every element to remove
        try:
            Lc.remove(rl)
        except ValueError:
            pass
    return Lc  # Return removed list

def Jaccard_Dis(s1,s2):
    # s1 and s2 are two sets.
    interS = list(set(s1) & set(s2))
    unionS =  list(set(s1) | set(s2))
    
    JaSim = float(len(interS))/len(unionS)
    JaSDis = 1 - JaSim
    return JaSDis
    
class CDBSCAN:
    # DBSCAN for discrite binary points
    def __init__(self, eps = 0.3, MinPts = 2, reuse_computed = 1):
        self.eps = eps
        self.MinPts = MinPts
        self.reuseC = reuse_computed
        
    def set_X(self, X):
        self.X = X
        self.Nsam = X.shape[0]  # Number of samples
        self.Ndim = X.shape[1]  # Number of datapoints
        
        ### Intialize Global Variables for the process
        # It assings to every point a cluster (-1 = Noise)
        self.cluster_P = -1 * np.ones((self.Nsam,1))  
        # Number of clusters so far. 
        self.K = -1;    
        self.samples_K = []  # Lists of lists. 
        # Each list inside constains the list of points of cluster i
        # We find new clusters when we find unassigned density points

        ################ OPTIMIZATION ###########
        # Number of elements in each sample. For obtaining the Jaccard distance later
        self.Nelem = self.X.getnnz(axis = 1)  
    
        if (self.reuseC == 1): # If we store precalculated distances
            self.Already_Calculated_Regions = [[-1]] * self.Nsam
            # When we calculate the region of a point, we also put here for the future.
            # If the samples is assigned to a cluster, then we remove it
        
        self.Nregion = 0;  # Number of times regions are calculated.
        # Just to have estimation of execution time.
        
    def regionQuery(self,i, mode = 1, only_same = 1):
        # This functions outputs the neighbourPoints of Point P given the threshold eps
        # With the variable mode we can choose the way the Region is calculated
        # With the variable only_same we can choose that we only consider neighbours
        # for the cluster, those points that havent been assigned to any cluster yet.
    
        eps = self.eps
        Nsam = self.Nsam
        
        if (self.reuseC == 1): # If we store precalculated distances
            ## If the region is precalculated !
            if (self.Already_Calculated_Regions[i][0] != -1):
                RegionSamples = self.Already_Calculated_Regions[i]
                print "Region Computation reused !!!"
                return RegionSamples
                
        # If we dont have precomputed ! We compute it :)

        self.Nregion += 1    # Regions calculated so far.
        if (self.Nregion % 50 == 0):   # Every 50 regions computed we print.
            print "Nregion: " + str(self.Nregion)
        
        ## Firt mode to calculate the Region
        if (mode == 0):  # Not efficient mode
            RegionSamples = []
            for j in range(Nsam):               # For every sample
                JaDis = self.JaDis_samples(i,j) # Get its metric
                if (JaDis <= eps):          # If the point j is in the region of i 
                    RegionSamples.append(j)
                    
        ### Efficient way to calculate the Jaccard distance
        # This function uses compresses sparse functions.
        elif (mode == 1):
            # Obtain the insertection by multiplicating the sample by all samples
            # and summing the results.
            inters = self.X[i,:].multiply(self.X[:,:]).sum(axis = 1)
            # Obtain the unions by means of the intersection and the number of samples.
            unions = float(self.Nelem[i]) + self.Nelem - inters.T
            # Compute JaDis_all to all samples.
            JaDis_all = 1 - inters /unions.T
            # Only get those ones closer than the threshols
            RegionSamples = np.where(JaDis_all <= self.eps )[0].tolist()
            
        if (self.reuseC == 1): # If we store precalculated distances
            self.Already_Calculated_Regions[i] = RegionSamples
        
        if (only_same == 0):  # If we just give all the neighbours
            return RegionSamples        
                
        #########################################################  
        ##### Optimization  ############################  
        ########################################################  
        # We could retrieve only the unassigned neighbors and the neighbours 
        # of the same cluster as the samples. This will be the real neighbours
        # to consider the point as a density point.
        cluster_P_Region = self.cluster_P[RegionSamples]
#        print cluster_P_Region
        non_assigned_indx = np.where(cluster_P_Region == -1)[0]
        
#        print non_assigned_indx
        non_assigned_points = np.array(RegionSamples)[non_assigned_indx].tolist()
        if (self.cluster_P[i] != -1):
            same_cluster_indx = np.where(cluster_P_Region == self.cluster_P[i])[0]
            same_cluster_points = np.array(RegionSamples)[same_cluster_indx].tolist()
            if(len(same_cluster_points) > 0 ):
                non_assigned_points.extend(same_cluster_points)

        return non_assigned_points
    

    def fit(self, X):  # Fit points X
        # This function obtaines the clusters for the compressed sparce matrix X

        self.set_X(X) # Initialize data structures
        
        for i in range(self.Nsam):  # For every sample
        # If point already assigned to cluster, we dont analyse it
          if (self.cluster_P[i] != -1):  
             continue;  # Get out of the loop and process next point
          
          # Get the samples that are inside the region of the point
          NeighborPts = self.regionQuery(i)  
          
          # If there are not enough data points
          if (len(NeighborPts) < self.MinPts):  
             self.cluster_P[i] = -1;  # Mark as noise for now
          
          else:  # Otherwise we start a Breath First Search of the nodes in the cluster
             self.K += 1;    # Increase the number of clusters.
             self.cluster_P[i] = self.K  # Set the point to the cluster
             
             if (self.reuseC == 1): # If we store precalculated distances
                 self.Already_Calculated_Regions[i] = [-1] # Remove its precalcualted neighbours
                 gc.collect()
             print "New cluster density point found: " + str(i)
             print "Obtaining cluster..."
             
             # Expand cluster !! Look for the nodes belonging to that cluster.
             self.expandCluster(i, NeighborPts)
             print "Cluster finished " + str(self.K) + " found: "+ \
                 str(len(self.samples_K[self.K])) + " samples."
        
        if (self.K == -1):
            print "No Clusters found" 
        
        ## Obtain the samples that outliers.
        self.K += 1
        outliers = np.where(self.cluster_P == -1)[0].tolist()
        self.samples_K.append(outliers)
    
    def print_clusters_sizes(self):
        # This function prints the sizes of all clusters
        size_clus = []
        for clus in self.samples_K:
            size_clus.append(len(clus))
        print "Sizes of the clusters"
        print size_clus 
        
    def JaDis_samples (self,i,j, mode = 1):
        # Jaccard Distance between two samples i and j in an unefficient way
        # self.X[i,:].todense()  # To convert to dense
        # Shouldnt be used in the final representation.
        if (mode == 0):
            set_i = np.where(self.X[i,:].toarray() == 1)[1] # Indexes of dimensions == 1
            set_j = np.where(self.X[j,:].toarray() == 1)[1]
            JaDis = Jaccard_Dis(set_i,set_j)  
        
        if (mode == 1):
            x1 = self.X[i,:].toarray()
            x2 = self.X[j,:].toarray()
            union = np.sum(np.logical_or(x1, x2))
            inter = x1.dot(x2.T)
            JaDis = 1- float(inter)/union
        return JaDis
            
    def expandCluster(self, i, NeighborPts):
        # Given a DensityPoint and its neighbours, this function
        # gets all the points that belong to the cluster.
        # NeighborPts is the initial set of points to check
       K_samples = []; # Indexes of the points that belong to the cluster
       K_samples.append(i)
       
       N_neighbors = len(NeighborPts) 
       # Number of total neighBours we have to look at.
       # If the neighbour has not been visited we add it to the cluster
       # If it is a density point, we also add its neighbours to this set
       j = 0;  # Index on the neighbors. It will iterate over new poins as well
       
       # Breath First Search of the nodes in the cluster  !!!     
       while (j < N_neighbors):
          Neigh_i = NeighborPts[j]
           # For every neighbour point (sample) to the density point
          
          # If a neighbour is in a cluster, it means that all its neighbours
          # are also already in a cluster... So no need to recalculate this.
          if (self.cluster_P[Neigh_i] == -1): # if P' neighbour does not belong to any clsutert 
             Neighb_ij = self.regionQuery(Neigh_i)  # Obtain the neighbour points 
             
             if (len(Neighb_ij) >= self.MinPts):  # If this point is also a density point
                # We add its Neighbours to the NeighborPts.
                # Only the neighbour that we were not considering already
                # We get the points in common between these new neighbors and the previous ones
                 InterNeigh = list(set(NeighborPts) & set(Neighb_ij))
                 NewNeigh = remove_from_list(Neighb_ij,InterNeigh)
                 if (len(NewNeigh) > 0): # If new neighbors found
                     NeighborPts.extend(NewNeigh)
#                 print NeighborPts
                 N_neighbors = len(NeighborPts)  # Recalculate number of Neighbours
            
             # if P' is not yet member of any cluster
             if ((self.cluster_P[Neigh_i]) == -1):  #
                 K_samples.append(Neigh_i)   # We assing it to the current cluster
                 self.cluster_P[Neigh_i] = self.K
                 # We remove the precalculated region !!
                 if (self.reuseC == 1): # If we store precalculated distances
                     self.Already_Calculated_Regions[Neigh_i] = [-1]
                     gc.collect()
          j+= 1; # Increase the index
        # Now we append all the samples of the cluster to the global structure
       self.samples_K.append(K_samples)

