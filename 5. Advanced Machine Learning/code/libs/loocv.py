import data_preprocessing as dp
import system_modules as sm
import numpy as np
import EM_libfunc as EMlf
import pandas as pd
from sklearn.model_selection import train_test_split


def load_data(filenames,path,conditions):
    """load all the .mat files in filenames.
    Args:
        filenames: list of .mat file names
        path: where files are locates
        conditions: list with classes to return

    Returns:
        X_subjects: list of subjects. Each subject has a list of classes. Each class is a
                    numpy array with dimensions [trials,samples,channels]
        Y_subjects: list of subject-labels. Each element subject-labels is a list of labels

    """
    X_subjects = []
    Y_subjects = []
    for subject in filenames:
        subject_data, subject_labels = dp.load_one_person_dataset(dataset_folder = path,
                                                              filename = subject)
        subject_selected_conditions = [subject_data[0],subject_data[2]]
        subject_selected_labels = [subject_labels[0],subject_labels[2]]
        # for condition in conditions:
        #     subject_selected_conditions.append(subject_data[condition])
        #     subject_selected_labels.append(subject_data[condition])

        X_subjects.append(subject_selected_conditions)
        Y_subjects.append(subject_selected_labels)

    return X_subjects,Y_subjects

def preprocess(X_subjects,Y_subjects, Nclasses, normalize = False):
    """Preprocessing and calculate the mean-trial of each subject
    Args:
        X_subjects: list of subjects. Each subject has a list of classes. Each class is a
                    numpy array with dimensions [trials,samples,channels]
        Y_subjects: list of subject-labels. Each element subject-labels is a list of labels
        Nclasses: number of classes
        normalize: ---
    Returns:
        X_averaged: List of subjects. Each subject is a len(classes) list. Each class contains
                    the averaged trial for that class.
    """
    X_subjects_averaged = []
    for i,subject in enumerate(X_subjects):
        labels = Y_subjects[i]

        channel_sel = range(70) #all channels
        # channel_sel = [20,21,22]
        max_trials = 500
        # THIS IS: channel selection, channel normalization, MAKE MODULUS 1?WHY?
        # we get a list: len=592. [0:295] class0 [296:591] class 1
        subject_trials, subject_labels = sm.create_data(subject, labels,
                    channel_sel = channel_sel, max_trials = max_trials,
                    rem_timePoint_ave = True, rem_features_ave = False,
                    normalize = False )

        average_subject = dp.get_average_from_train(Nclasses, subject_trials, subject_labels)

        # ##### Separate in train and validation
        # r_seed = np.abs(int(100 * np.random.randn()))
        #
        # X_train, X_test, y_train, y_test = train_test_split(subject_trials,
        #             subject_labels, test_size=0.50, random_state = r_seed, stratify = subject_labels)
        #
        # if(normalize == True):
        #     X_train = dp.normalize_trialList(X_train)
        #     X_test = dp.normalize_trialList(X_test)
        #
        # X_subjects_train.append(X_train)
        # X_subjects_test.append(X_test)
        # Y_subjects_train.append(y_train)
        # Y_subjects_test.append(y_test)
        X_subjects_averaged.append(average_subject)

    return X_subjects_averaged

def crossvalidation_EM(X_subjects, Klist, saveInfo = True, path = 'file.csv'):
    """LOO crossvalidation of EM for number of clusters and each class
    Args:
        X_subjects: list of subjects with from from one class
        Klist: list of number of clusters to crossvalidate
    Returns:
        likelihood_list_train:  list of likelihoods over training for each K
        likelihood_list_test: list of likelihoods over testing for each K
    """
    Nclasses = len(X_subjects[0])
    Nsubjects = len(X_subjects)
    nanVector = np.full(len(Klist),np.nan)
    dic = {'Klist':Klist,'train':nanVector,'test':nanVector}; df = pd.DataFrame(dic)
    for i,K in enumerate(Klist):
        average_likelihood_train = 0
        average_likelihood_test = 0
        for idx in range(Nsubjects):
            validation = X_subjects[idx]
            train = cutout(X_subjects,idx)
            Ninit = 1; maxIt = 100;
            Ks_params = sm.get_clusters_labels_EM(1, X_train = train, y_train = [0] * len(train),
                                     Ninit = Ninit, K = K, T = maxIt ,verbose = 0)
            #Evaluate likelihood
            likelihood_train = EMlf.get_EM_Incomloglike_log(Ks_params[0][1],Ks_params[0][0],X = concatenateTrials(train))
            likelihood_test = EMlf.get_EM_Incomloglike_log(Ks_params[0][1],Ks_params[0][0],X = validation)
            average_likelihood_train += likelihood_train
            average_likelihood_test += likelihood_test
        average_likelihood_train = average_likelihood_train/(float(Nsubjects)*(Nsubjects-1))
        average_likelihood_test = average_likelihood_test/float(Nsubjects)
        df.set_value(i, 'test', average_likelihood_test)
        df.set_value(i, 'train', average_likelihood_train)
        if saveInfo==True:
            df.to_csv(path)
    return df

def crossvalidation_HMM(X_subjects, Klist, saveInfo = True, path = 'file.csv'):
    """LOO crossvalidation of HMM for number of clusters and each class
    Args:
        X_subjects: list of subjects with from from one class
        Klist: list of number of clusters to crossvalidate
    Returns:
        likelihood_list_train:  list of likelihoods over training for each K
        likelihood_list_test: list of likelihoods over testing for each K
    """
    Nclasses = len(X_subjects[0])
    Nsubjects = len(X_subjects)
    nanVector = np.full(len(Klist),np.nan)
    dic = {'Klist':Klist,'train':nanVector,'test':nanVector}; df = pd.DataFrame(dic)
    for i,K in enumerate(Klist):
        average_likelihood_train = 0
        average_likelihood_test = 0
        for idx in range(Nsubjects):
            validation = X_subjects[idx]
            x_train = cutout(X_subjects,idx)
            y_train = [0] * len(x_train)
            Nit = 2; I = K; maxIt = 100; Ks_params_HMM = None
            Is_params = sm.get_clusters_labels_HMM(1,X_train = x_train, y_train = y_train,
                            Ks_params = Ks_params_HMM, Nit =Nit, I  =  I, R  = maxIt, verbose = 1)
            #Evaluate likelihood
            like_train = dp.get_likelihoods_HMM(x_train,Is_params); like_train = sum(like_train)
            like_test = dp.get_likelihoods_HMM([validation],Is_params); like_test = sum(like_test)
            average_likelihood_train += like_train[0]
            average_likelihood_test += like_test[0]
        average_likelihood_train = average_likelihood_train/(float(Nsubjects)*(Nsubjects-1))
        average_likelihood_test = average_likelihood_test/float(Nsubjects)
        df.set_value(i, 'test', average_likelihood_test)
        df.set_value(i, 'train', average_likelihood_train)
        if saveInfo==True:
            df.to_csv(path)
    return df

def concatenateTrials(X):
    """Concatenate trials in X and returns a numpy X_concatenated"""
    X_concatenated = X[1][:]
    for trial in X[1:]:
        X_concatenated = np.concatenate((X_concatenated,trial))
    return X_concatenated

def divideClass(subjects):
    """Create a list per class
    Args:
        subjects: list of subjects. Each subject has N classes
    return:
        subjects_classes: return a list for each class
    """
    subject_classes = []
    Nclasses = len(subjects[0])
    for i in range(Nclasses):
        subject_classes.append([])
    for subject in subjects:
        for j,condition in enumerate(subject):
            subject_classes[j].append(condition)
    return subject_classes

def cutout(seq, idx):
    """ Remove element at 'idx' from 'seq' """
    return seq[:idx] + seq[idx + 1:]

def get_average_trials(X_train,Y_train,Nclasses):
    """For each subjects all the trials are averaged"""
    subjects_mean_trial = []
    for i in range(Nclasses):
        class1 = []
        for j,subject in enumerate(X_train):
            X_data_ave = dp.get_average_from_train(Nclasses, subject, Y_train[j])
            class1.append(X_data_ave[i])
        subjects_mean_trial.append(class1)
    return subjects_mean_trial
