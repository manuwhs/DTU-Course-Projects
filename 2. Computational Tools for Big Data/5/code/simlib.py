#### SIMILARITY LIBRARY #####

import numpy as np
from utilfunc import *
from textprocesslib import *
# Define similarity values between the query and the corpus.
def get_similarities(query,database):
     query_words = doc_preprocess(query,1)  # Preprocess query
     query_words_tfidf = tf_idf_doc(query_words, database) # Get tfidf of query
     
     for i in range (len(query_words_tfidf[0])):
         print str(query_words_tfidf[1][i]) + "\t\t" + query_words_tfidf[0][i]
     print " "
     N_docs_database = len(database)
     
     eu_dists = np.zeros([N_docs_database,1])
     cos_sims = np.zeros([N_docs_database,1])
     n_commun = np.zeros([N_docs_database,1])
     tfidf_sum = np.zeros([N_docs_database,1])
     
     for i in range (N_docs_database):
         # Common words
         (vq ,vd) = get_common_tfidf_vectors(query_words_tfidf,database[i])
         # Non-common words
         (vnq ,vnd) = get_noncommon_tfidf_vectors(query_words_tfidf,database[i])
         
         # If no words in commun
         if(vq.size == 0):
              eu_dists[i] = -1
              cos_sims[i] = -1
              n_commun[i] = vq.size
         else:
             eu_dists[i] = get_euc_dist(vq,vd,vnq,vnd)
             cos_sims[i] = get_cos_sim(vq,vd)
             n_commun[i] = vq.size
             tfidf_sum[i] = np.sum(vd)
             
     similarities = (eu_dists,cos_sims, n_commun, tfidf_sum)
     # Similatiries obtained are:
     # - The euclidean distance
     # - The cosine similarity
     # - The number of words in commun. 
     #      (Important coz, they can have very high values in the eu y cos
     #      but just because they have very few words in commun)
     # - tfidf_sum: Sum of the tdifd of the values
     
     return similarities

def rank_documents(Similarity, N_top, type_simi):
    # RANK WITH EVERY INDEPENDENT SIMILARITY MEASURE.
    n_doc = len(Similarity)
    # 1) Euclidean distance: The closer, the better, so the smaller the value, the better
    if (type_simi == "euclidean"):
        # For this similarity, the vectors that did not have words in common have
        # similarity -1, so we cannot order them directly like this. 
        for i in range(len(Similarity)):
            if (Similarity[i] == -1):
                Similarity[i] = 10000000000;
        indexes = np.array(range(0,n_doc))
        together = zip(Similarity, indexes)  # Zip both arrays together
        sorted_together =  sorted(together) # Sorted Increasing order
        ordered_indexes = [x[1] for x in sorted_together]
        # Since ordered_indexes contains the indexes of the best documents 
        # (obtained from decreasing similarity)
        
    # 2) Cosine distance: The biger the value, the better 
    if (type_simi == "cosine"):
        # Wrong vectors have similarity -1 so no problem, the bigger the better.
        indexes = np.array(range(0,n_doc))
        together = zip(Similarity, indexes)  # Zip both arrays together
        sorted_together =  sorted(together, reverse=True) # Sorted Decreasing order
        ordered_indexes = [x[1] for x in sorted_together]
        #print ordered_indexes
    
    # 3) Use the number of words in common:
    if (type_simi == "common"):
        indexes = np.array(range(0,n_doc))
        together = zip(Similarity, indexes)  # Zip both arrays together
        sorted_together =  sorted(together, reverse=True) # Sorted Decreasing order
        ordered_indexes = [x[1] for x in sorted_together]
        
    # ) Use the sum of tfidf values:
    if (type_simi == "tdidf_sum"):
        indexes = np.array(range(0,n_doc))
        together = zip(Similarity, indexes)  # Zip both arrays together
        sorted_together =  sorted(together, reverse=True) # Sorted Decreasing order
        ordered_indexes = [x[1] for x in sorted_together]
        
    top_indexes = ordered_indexes[0:N_top]
    return top_indexes

def get_combined_rank(Ranks, N_top, type_comb):
# Ok... we have several similarity measures, now we have to output the ranking
# based on those similarities. How do we combine them to get the best ranking ?
    n_rank,n_doc = np.shape(Ranks)
    
    eu_w = 0.5; # Weight of the euclidean distance importance
    cos_w = 0.5; # Weight of the cosine distance importance
    n_w = 1;   # Weight of the nomber of common words importance
    
    total_similarity  = np.zeros([n_doc,1])
    
    for i in range(n_doc):
        total_similarity [Ranks[i]] += eu_w * Ranks[0][i] 
        total_similarity [Ranks[i]] += cos_w * Ranks[1][i] 
        total_similarity [Ranks[i]] -= n_w * Ranks[2][i] 
    
    # Get the top indexes rank from the combined similarity methods:
    indexes = np.array(range(0,n_doc))
    together = zip(total_similarity, indexes)  # Zip both arrays together
    sorted_together =  sorted(together) # Sorted Increasing order
    ordered_indexes = [x[1] for x in sorted_together]
    top_indexes = ordered_indexes[0:N_top]
    return top_indexes
    
def get_euc_dist(vq,vd,vnq,vnd):
    distance = np.sqrt(((vq - vd)*(vq - vd)).sum(axis=0))

    
    return distance
    
def get_cos_sim(vq,vd):
    modules = (np.sqrt((vq*vq).sum(axis=0)) * np.sqrt((vd*vd).sum(axis=0)))
    
    if (modules == 0):
        print "Error, cosine similarity. One vector is all 0s."
        return -1
    cos_sim = np.dot(vq.transpose(),vd)
    cos_sim = cos_sim / modules
    
    return cos_sim  
