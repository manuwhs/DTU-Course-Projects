
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os as os
import matplotlib.colors as ColCon
from scipy import spatial
import datetime as dt

w = 10  # Width of the images
h = 6   # Height of the images

# Define the empty dataframe structure
keys = ['Open', 'High', 'Low', 'Close', 'Volume']
empty_df= pd.DataFrame(None,columns = keys )

keys_col = ['Symbol','Type','Size','TimeOpen','PriceOpen', 'Comision','CurrentPrice','Profit']
empty_coliseum = pd.DataFrame(None,columns = keys_col )


# Dictionary between period names and value
periods = [1,5,15,30,60,240,1440,10080,43200, 43200*12]
periods_names = ["M1","M5","M15","M30","H1","H4","D1","W1","W4","Y1"]
period_dic = dict(zip(periods,periods_names))
names_dic = dict(zip(periods_names, periods))

def fnp(ds):
    # This function takes some numpy element or list and transforms it
    # into a valid numpy array for us.
    # It works for lists arrays [1,2,3,5], lists matrix [[1,3],[2,5]]
    # Vectors will be column vectors
    
    # Working with lists
    if (type(ds).__name__ == "list"):
        # If the type is a list 
        if (len(ds) == 0):  # If we are given an empty list 
            ds = np.array(ds).reshape(1,0)
            return ds
        
        ## If what is inside is a vector
        if(np.array(ds[0]).size > 1):
            ds = np.array(ds)
            ds = ds.reshape(np.array(ds).size/np.array(ds[0]).size,np.array(ds[0]).size)
            return ds
        else:
            ds = np.array(ds)
#            print ds.shape
            ds = ds.reshape(ds.size,1) # Return column vector
    
    # Working with nparrays
    elif (type(ds).__name__ == 'numpy.ndarray' or type(ds).__name__ == "ndarray"):

        if (len(ds.shape) == 1): # Not in matrix but in vector form 
            ds = ds.reshape(ds.size,1)
            
        elif(ds.shape[0] == 1):
            # If it is a row vector instead of a column vector.
            # We transforme it to a column vector
            ds = ds.reshape(ds.size,1)
            
    elif (type(ds).__name__ == 'DatetimeIndex'):
        ds = pd.to_datetime(ds)
        ds = np.array(ds).reshape(len(ds),1) 
        
    return ds
    
def convert_to_matrix (lista, max_size = -1):
    # Converts a list of lists with different lengths into a matrix 
    # filling with -1s the empty spaces 

    Nlist = len(lista)
    
    listas_lengths = []
    
    if (max_size == -1):
        for i in range (Nlist):
            listas_lengths.append(lista[i].size)
        
        lmax = np.max(listas_lengths)
    else:
        lmax = max_size 
        
    matrix = -1 * np.ones((Nlist,lmax))
    
    for i in range (Nlist):
        if (lista[i].size > lmax):
            matrix[i,:lista[i].size] = lista[i][:lmax].flatten()
        else:
            matrix[i,:lista[i].size] = lista[i].flatten()
    
    return matrix
    
def get_dates(dates_list):
    # Gets only the date from a timestapm. For a list
    only_day = []
    for date in dates_list:
        only_day.append(date.date())
    return np.array(only_day)
    
def str_to_datetime(dateStr):
    # This function converts a str with format YYYY-MM-DD HH:MM:SS to datetime
    dates_datetime = []
    for ds in dateStr:
        dsplited = ds.split(" ")
        date_s = dsplited[0].split("-") # Date
        
        if (len(dsplited) > 1):  # Somo files have hours, others does not
            hour_s = dsplited[1].split(":")  # Hour 
            datetim = dt.datetime(int(date_s[0]), int(date_s[1]), int(date_s[2]),int(hour_s[0]), int(hour_s[1]))
        else:
            datetim = dt.datetime(int(date_s[0]), int(date_s[1]), int(date_s[2]))
            
        dates_datetime.append(datetim)
    return dates_datetime
    
def load_dataset(file_dir = "./dataprices.csv"):
    # Reads the dataprices
    data = pd.read_csv(file_dir, sep = ',') # header = None, names = None  dtype = {'phone':int}
    Nsamples, Ndim = data.shape   # Get the number of bits and attr

    return data

def sort_and_get_order (x, reverse = True ):
    # Sorts x in increasing order and also returns the ordered index
    x = x.flatten()  # Just in case we are given a matrix vector.
    order = range(len(x))
    
    if (reverse == True):
        x = -x
        
    x_ordered, order = zip(*sorted(zip(x, order)))
    
    if (reverse == True):
        x_ordered = -np.array(x_ordered)
        
    return np.array(x_ordered), np.array(order)
    
def create_folder_if_needed (folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
    

def windowSample (sequence, L):
    """ Transform a sequence of data into a Machine Learning algorithm,
    it transforms the sequence into X and Y being """
    
    sequence = np.array(sequence).flatten()
    Ns = sequence.size
    
    X = np.zeros((Ns - (L +1), L ))
    Y = np.zeros((Ns - (L +1),1) )
    for i in range (Ns - (L +1)):
        X[i,:] = sequence[i:i+L]
        Y[i] = sequence[i+L]
    # We cannot give the output of the first L - 1 sequences (incomplete input)
    return X, Y

def simmilarity(patterns,query,algo):
    # This funciton computes the similarity measure of every pattern (time series)
    # with the given query signal and outputs a list of with the most similar and their measure.

    Npa,Ndim = patterns.shape
    sims = []
    if (algo == "Correlation"):
        for i in range(Npa):
            sim =  np.corrcoef(patterns[i],query)[1,0]
            sims.append(sim)
        sims = np.array(sims)
        sims_ored, sims_or = sort_and_get_order (sims, reverse = True )
        
    if (algo == "Distance"):
        sims = spatial.distance.cdist(patterns,np.matrix(query),'euclidean')
        sims = np.array(sims)
        sims_ored, sims_or = sort_and_get_order (sims, reverse = False )
    return sims_ored, sims_or

def get_Elliot_Trends (yt, Nmin = 4, Noise = -1):
    
    Nsamples, Nsec = yt.shape
    if (Nsec != 1):
        print "Deberia haber solo una senal temporal"
        return -1;
        
#    yt = yt.ravel()
    
#    yt = np.array(yt.tolist()[0])
    
    print yt.shape
    trends_list = []   # List of the trends
    
    support_t = 0   # Support index
    trend_ini = 0   # Trend start index

    support = yt[support_t]  # If support is broken then we dont have trend
    

    """ UPPING TRENDS """    
    for i in range (1,Nsamples-1):
        if (Noise == -1):
            tol = support/200
            
        #### Upper trends
        if (yt[i] > support- tol): # If if is not lower that the last min
            if (yt[i +1 ] < yt[i] - tol):  # If it went down, we have a new support
                support_t = i
                support = yt[support_t]
            
        else:   # Trend broken
            
            if ((i -1 - trend_ini) > Nmin): # Minimum number of samples of the trend
                trends_list.append([trend_ini, i -1])  # Store the trend
            
            # Start over
            trend_ini = i
            support_t = i
            support = yt[support_t]
    
    """ Lowing TRENDS """  
    
    for i in range (1,Nsamples-1):
        if (Noise == -1):
            tol = support/200
            
        #### Upper trends
        if (yt[i] < support + tol): # If if is not lower that the last min
            if (yt[i + 1] > yt[i] + tol):  # If it went up, we have a new support
                support_t = i
                support = yt[support_t]
            
        else:   # Trend broken
            
            if ((i - trend_ini) > Nmin): # Minimum number of samples of the trend
                trends_list.append([trend_ini, i -1])  # Store the trend
            
            # Start over
            trend_ini = i
            support_t = i
            support = yt[support_t]
    return trends_list
        
def plot_trends(yt, trends_list):
    
    for trend in trends_list:
        plt.plot(trend, yt[trend], lw = 5)
        

def support_detection(sequence, L):
    # This fuction get the support of the last L signals
    Nsamples, Nsec = sequence.shape
    
    sequence_view = sequence[-L:]
    index_min = np.argmin(sequence_view)
    
    return index_min + Nsamples - L 