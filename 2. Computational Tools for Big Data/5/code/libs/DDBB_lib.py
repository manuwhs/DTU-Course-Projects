# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 23:55:37 2016

@author: montoya
"""
#

import pandas as pd
import numpy as np
import urllib2
import datetime as dt
import os as os 

def create_folder_if_needed (folder):
    if not os.path.exists(folder):
        os.makedirs(folder)

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
    
def save_to_csv(symbol,dataCSV, file_dir = "./storage/"):
    create_folder_if_needed(file_dir)
    whole_path =  file_dir + symbol + ".csv"
    dataCSV.to_csv(whole_path, sep=',')

def load_csv_timeData(symbol, file_dir = "./storage/"):

    whole_path = file_dir + symbol + ".csv"
    try:
        dataCSV = pd.read_csv(whole_path,
                          sep = ',', index_col = 0, dtype = {"Date":dt.datetime})
    
        dataCSV.index = str_to_datetime (dataCSV.index.tolist())
        
    except IOError:
        error_msg = "File does not exist: " + whole_path 
        print error_msg
    except:
        print "Unexpected error in file: " + whole_path
    # We transform the index to the real ones
    return dataCSV
    

def get_data_yahoo(symbol = "AAPL", precision = "1mo", 
                   start_date = "01-12-2011", end_date = "01-12-2015"):

    # data1 = dt.datetime.fromtimestamp(1284101485)
    sdate = dt.datetime.strptime(start_date, "%d-%m-%Y")
    edate = dt.datetime.strptime(end_date, "%d-%m-%Y")
    
#    sdate_ts = int(get_timeStamp(sdate))
#    edae_ts = int(get_timeStamp(edate))

#    url_root = "https://finance.yahoo.com/quote/"
#    url_root += symbol
#    url_root += "/history?"
#    url_root += "period1=" + str(sdate_ts)
#    url_root += "&period2=" + str(edate_ts)
#    url_root += "&interval=" + precision
#    url_root += "&filter=history&frequency=" + precision
 
    url_root = "http://chart.finance.yahoo.com/table.csv?"
    url_root += "s=" + symbol
    url_root += "&a=" +str(sdate.day)+ "&b=" +str(sdate.month)+ "&c=" +str(sdate.year)
    url_root += "&d=" +str(edate.day)+ "&e=" +str(edate.month)+"&f="+str(edate.year)
    url_root += "&g=" + "m"
    url_root += "&ignore=.csv"
    
#    print url_root
    response = urllib2.urlopen(url_root)
    data = response.read().split('\n')
    nlines = len(data)
    for i in range(nlines):
        data[i] = data[i].split(",")
        
#    print data[0:4]
    df = pd.DataFrame(data)
    df.columns = df.ix[0]  #['Date','Open', 'High', 'Low', 'Close', 'Volume', "Adj Close"]
#    print df.columns
 
    
    ### REMOVE FIRST ROW (Headers) 
    df.drop(0, inplace = True)
    ### REMOVE LAST ROW (Nones)
#    print len(df) - 1
#    print df.ix[len(df) - 1]
    df.drop(df.index.values[len(df) - 1], inplace = True)
    ### CONEVERT DATES TO TIMESTAMPS (Nones)
#    print df.Date
    df.index = str_to_datetime(df.Date)
    
    del df['Date']
    ## 
    # We have to
    return df


def download_and_add(list_symbols, sdate = "01-01-1996",
                     edate = "01-01-2016", fir_dir =  "./storage"):
    
    for symbol in list_symbols:
        data_Symbol = get_data_yahoo (symbol = symbol,  precision = "1mo", 
                  start_date = sdate, end_date = edate)
                  
        save_to_csv(symbol = symbol, 
                    dataCSV = data_Symbol)
