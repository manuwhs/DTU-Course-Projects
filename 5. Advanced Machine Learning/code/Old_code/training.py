# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 01:31:58 2015

@author: montoya
"""

import numpy as np
from sklearn import cross_validation
import paramClasses as paC

from time import time    # For the random seed

def train (self):
 
    # Adapt the labels so that they are correct (-1 or 0 and transform multivariate if needed)
    self.Ytrain = paC.adapt_labels(self.Ytrain, mode = self.outMode)
    self.Yval = paC.adapt_labels(self.Yval, mode = self.outMode )
    if (self.Xtest != []):     # If there is a test dataset.
        self.Ytest = paC.adapt_labels(self.Ytest, mode = self.outMode )
    
    
    # Values that will be stored in Exec_list for later processing
    self.TrError = np.zeros((self.Nruns,1))
    self.ValError = np.zeros((self.Nruns,1))
    self.TstError = np.zeros((self.Nruns,1))
    
    for r in range (self.Nruns):
        self.train_CV(r = r)
        
#        print self.TrError[r], self.ValError[r] , self.TstError[r]

def train_once (self):    
#    print "train_once"
#    print self.D
#    print self
    
    self.init_Weights()                  # Initialize

    # Check the training algorithm and pass it with its parameters.
    # D is the dehenfasis vector, distribution of samples probabilities.
    if (self.trainingAlg.trAlg == "ELM"):
        self.ELM_train(self.trainingAlg.param)
        
    if (self.trainingAlg.trAlg == "BP"):
        self.BP_train(self.trainingAlg.param) 
        
    if (self.trainingAlg.trAlg == "BMBP"):
        self.BMBP_train(self.trainingAlg.param) 
        
    if (self.trainingAlg.trAlg == "ELMT"):
        self.ELMT_train(self.trainingAlg.param)     

    if (self.trainingAlg.trAlg == "LDAT"):
        self.LDAT_train(self.trainingAlg.param)     

    if (self.trainingAlg.trAlg == "LDA"):
        self.LDA_train(self.trainingAlg.param)    
            
def train_CV (self, r):
    # Trains the learner CV times using cross validation

    total_Xtrain = self.Xtrain
    total_Ytrain = self.Ytrain
    ## Get the random seed and use it
    if (self.InitRandomSeed == -1): # If no seed is specified
        self.RandomSeed[r] = int((time()%1 * 100000))
        np.random.seed(self.RandomSeed[r])
    else:
        self.RandomSeed[r] = self.InitRandomSeed
        np.random.seed(self.RandomSeed[r])
        
    TrError = 0;
    ValError = 0;
    TstError = 0;
    
#    print "train_CV"
#    print self.CV
#    print self.D
    
    if (self.CV == 1):  
        # If the validation is performed with just the training set
        # Then the validation set is the original self.Xval. 
        """ Why you may ask ?? """ 
        # In other aggregate solutions, like Boosting, the CV is done
        # over the whole structure, not layer by layer. In this cases,
        # the CV of the SLFN will be 1 always and its the Boosting "train"
        # the one in charge for changing the Validation set and training set.

        self.train_once()
        
        TrError += self.score(self.Xtrain, self.Ytrain)
        ValError += self.score(self.Xval, self.Yval)
        if (self.Xtest != []):     # If there is a test dataset.
            TstError += self.score(self.Xtest, self.Ytest)
    
    if (self.CV > 1):

        stkfold = cross_validation.StratifiedKFold(total_Ytrain.ravel(), n_folds = self.CV)
        for train_index, val_index in stkfold:
#                print train_index
            self.set_Train(total_Xtrain[train_index],total_Ytrain[train_index])
            self.set_Val(total_Xtrain[val_index],total_Ytrain[val_index])
            self.train_once()
        
            TrError += self.score(self.Xtrain, self.Ytrain)
            ValError += self.score(self.Xval, self.Yval)
            if (self.Xtest != []):     # If there is a test dataset.
                TstError += self.score(self.Xtest, self.Ytest)
            
    self.TrError[r] = TrError / self.CV
    self.ValError[r] = ValError / self.CV
    self.TstError[r] = TstError / self.CV
    
    self.Xtrain = total_Xtrain   # Restore the original Xtrain
    self.Ytrain = total_Ytrain 
    
    