#########################################################3
############### DDBB LIBRARY  ##############################
##########################################################
## Library with function to manage the databases

import pandas as pd
import numpy as np
import urllib2
import datetime as dt
import utilities_lib as ul
import get_data_lib as gdl

def save_to_csv(symbol,dataCSV, file_dir = "./storage/"):
    ul.create_folder_if_needed(file_dir)
    whole_path =  file_dir + symbol + ".csv"
    dataCSV.to_csv(whole_path, sep=',')

def load_csv_timeData(symbol, file_dir = "./storage/"):

    whole_path = file_dir + symbol + ".csv"
    try:
        dataCSV = pd.read_csv(whole_path,
                          sep = ',', index_col = 0, dtype = {"Date":dt.datetime})
    
        dataCSV.index = ul.str_to_datetime (dataCSV.index.tolist())
        
    except IOError:
        error_msg = "File does not exist: " + whole_path 
        print error_msg
    except:
        print "Unexpected error in file: " + whole_path
    # We transform the index to the real ones
    return dataCSV

def load_dataset(file_dir = "./dataprices.csv"):
    # Reads the dataprices
    data = pd.read_csv(file_dir, sep = ',') # header = None, names = None  dtype = {'phone':int}
    Nsamples, Ndim = data.shape   # Get the number of bits and attr

    return data

def download_and_add(list_symbols, sdate = "01-01-1996",
                     edate = "01-01-2016", fir_dir =  "./storage"):
    
    for symbol in list_symbols:
        data_Symbol = gdl.get_data_yahoo (symbol = symbol,  precision = "1mo", 
                  start_date = sdate, end_date = edate)
                  
        save_to_csv(symbol = symbol, 
                    dataCSV = data_Symbol)
