
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
import scipy.io
from sklearn import preprocessing

import Watson_distribution as Wad
import Watson_sampling as Was
import Watson_estimators as Wae
import general_func as gf

import data_preprocessing as dp
import system_modules as sm

import loocv as art
plt.close("all")


############################ Load the dataset ! ##############################
filenames =  ['face_scrambling_spm8proc_sub01.mat','face_scrambling_spm8proc_sub02.mat',
            'face_scrambling_spm8proc_sub03.mat','face_scrambling_spm8proc_sub04.mat',
            'face_scrambling_spm8proc_sub05.mat','face_scrambling_spm8proc_sub06.mat',
            'face_scrambling_spm8proc_sub07.mat','face_scrambling_spm8proc_sub08.mat',
            'face_scrambling_spm8proc_sub09.mat','face_scrambling_spm8proc_sub10.mat',
            'face_scrambling_spm8proc_sub11.mat','face_scrambling_spm8proc_sub12.mat',
            'face_scrambling_spm8proc_sub13.mat','face_scrambling_spm8proc_sub14.mat',
            'face_scrambling_spm8proc_sub15.mat','face_scrambling_spm8proc_sub16.mat']
path = "./dataset/"
filename = [filenames[1]] # if passing a single subject
classes = [0,2]; #0:famous face, 1:unkown face, 2: scramble face
X_subjects, Y_subjects = art.load_data(filenames[:], path, classes)
Nclasses = len(classes)

########################### Preprocessing ###################################
X_subjects_averaged = art.preprocess(X_subjects,Y_subjects, Nclasses, normalize = False)
X = art.divideClass(X_subjects_averaged)

################################## CrossValidating ##############################
CVEM = 0; CVHMM = 1
Klist = [34]
if CVEM ==1:
    df_class0_EM = art.crossvalidation_EM(X[0], Klist, saveInfo = True, path='visualize/masEMclass0.csv')
    # df_class1_EM = art.crossvalidation_EM(X[1], Klist, saveInfo = True, path='visualize/EMclass1.csv')
if CVHMM ==1:
    df_class0 = art.crossvalidation_HMM(X[0], Klist, saveInfo = True, path='visualize/masHMMclass0.csv')
    # df_class1 = art.crossvalidation_HMM(X[1], Klist, saveInfo = True, path='visualize/HMMclass1.csv')

import IPython
IPython.embed()
