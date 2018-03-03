import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import numpy as np
from numpy import loadtxt
import time

import graph_lib as gr


def separate_days (price, date):
    # This function separates the unidimensional price data into a list of days.
    # date has the day.minutes
    
    date = np.array(date,dtype = int)  # We only keep the indicator of the day
    date_days = date/1000000;  # Vector of days 
    hour_days = date%1000000;  # Vector of time of the days
    
    days = np.unique(date_days)   # Different days
    Ndays = days.size
    
    Price_list_days = []
    Hours_list_days = []
    
    for i in range (Ndays):  # Get the price y hora of every day
        index_day = np.where(date_days == [days[i]])
        
        prices_day = price[index_day]
        hour_day = hour_days[index_day]
        
        Price_list_days.append(prices_day)
        Hours_list_days.append(hour_day)
        
    return Price_list_days, Hours_list_days, days

def time_normalizer (hour_data, prices, time_span_sec):
    # This function normalizes the sample vector given and creates a list
    # where every element contains the mean price of the timespan. 
    # If there is no time for a Span, we give it the price of the previous

    # Hours = HHMMSS 
    """ We can use the histogram function"""
    # Transform the hours into seconds
    seconds_data = (hour_data/10000)*3600 + ((hour_data/100)%100)*60 + hour_data%100

    price_day = []
    
    total_s = 24*60*60
    N_bins = total_s  /time_span_sec
#    print seconds_data[:100]
    for i in range (N_bins):
        bin_prices_indx = np.where(seconds_data < (i+1)*time_span_sec)
        bin_prices = prices[bin_prices_indx]
        
        Npr = bin_prices_indx[0].size
        seconds_data = seconds_data[Npr:]  # Remove the previous 
        prices = prices[Npr:]  # Remove the previous 
                                                    # assuming prices are ordered.
#        print bin_prices
        
        if (bin_prices.size == 0): # If we have no data for the slice
            if (i == 0):  # If it is the first slice
                price_day.append(-1)
            else:
                price_day.append(price_day[-1])
        else:
            price_day.append(np.mean(bin_prices))
    
    return price_day

def transform_time(time_formated):
    # This function accepts time in the format 2016-01-12 09:03:00
    # And converts it into the format [days] [HHMMSS]
    # Remove 
    
    data_normalized = []
    for time_i in time_formated:
        time_i = str(time_i)
#        print time_i
        time_i = time_i[0:19]
        time_i = time_i.replace("-", "")
        time_i = time_i.replace(" ", "")
        time_i = time_i.replace(":", "")
        time_i = time_i.replace("T", "")
#        print time_i
        data_normalized.append(int(time_i))
        
    return data_normalized 
        
def remove_list_indxs(lista, indx_list):
    # Removes the set of indexes from a list
    removeset = set(indx_list)
    newlist = [v for i, v in enumerate(lista) if i not in removeset]
    
    return newlist

def get_close_open_diff (prices_days):
    # Compute the difference between openning and closing prices
    # prices_days[days][prices]
    shape = prices_days.shape
    
    diff = []
    
    for i in range (shape[0] - 1):
        diff.append(prices_days[i+1][0] - prices_days[i][-1])
    
    diff = np.array(diff)
    return diff
    
def get_open_close_diff (prices_days):
    # Compute the difference between openning and closing prices of the same day
    # prices_days[days][prices]
    shape = prices_days.shape
    
    diff = []
    
    for i in range (shape[0]):
        diff.append(prices_days[i][-1] - prices_days[i][0])
    
    diff = np.array(diff)
    return diff
    
#==============================================================================
#     hours = hour_data/10000
#     mins = (hour_data/100)%100
#     
#     Allhours = np.unique(hours)
#     Nhours = Allhours.size
#     
#     N_samples_by_hour = 60/time_span_min
#     
#     price_day = []
#     for h in range (Nhours):
#         prices_indx_hour = np.where(hours == Allhours[h])  # Subselect the hour index
#         princes_hour = prices[prices_indx]
#         
#         for m in range (N_samples_by_hour):  # Subselect the spans (histogram)
#             prices_indx =  
#             price_span = np.where(prices[prices_indx] < time_span_min*
#             price_day.append()
#==============================================================================
        
        
        
def patternFinder():
    '''
    The goal of patternFinder is to begin collection of %change patterns
    in the tick data. From there, we also collect the short-term outcome
    of this pattern. Later on, the length of the pattern, how far out we
    look to compare to, and the length of the compared range be changed,
    and even THAT can be machine learned to find the best of all 3 by
    comparing success rates.
    '''
    
    #Simple Average
    avgLine = ((bid+ask)/2)
    
    #This finds the length of the total array for us
    x = len(avgLine)-30
    #This will be our starting point, allowing us to compare to the
    #past 10 % changes. 
    y = 11
    # where we are in a trade. #
    # can be none, buy,
    currentStance = 'none'
    while y < x:
        
        p1 = percentChange(avgLine[y-10], avgLine[y-9])
        p2 = percentChange(avgLine[y-10], avgLine[y-8])
        p3 = percentChange(avgLine[y-10], avgLine[y-7])
        p4 = percentChange(avgLine[y-10], avgLine[y-6])
        p5 = percentChange(avgLine[y-10], avgLine[y-5])
        p6 = percentChange(avgLine[y-10], avgLine[y-4])
        p7 = percentChange(avgLine[y-10], avgLine[y-3])
        p8 = percentChange(avgLine[y-10], avgLine[y-2])
        p9 = percentChange(avgLine[y-10], avgLine[y-1])
        p10= percentChange(avgLine[y-10], avgLine[y])

        outcomeRange = avgLine[y+20:y+30]
        currentPoint = avgLine[y]

        #function to account for the average of the items in the array
        print reduce(lambda x, y: x + y, outcomeRange) / len(outcomeRange)

        
        print currentPoint
        print '_______'
        print p1, p2, p3, p4, p5, p6, p7, p8, p9, p10
        time.sleep(55)
        
        y+=1
        


def graphRawFX():
    
    fig=plt.figure(figsize=(10,7))
    
    ax1 = plt.subplot2grid((40,40), (0,0), rowspan=40, colspan=40)
    ax1.plot(date,bid)
    ax1.plot(date,ask)
    
    ax1.plot(date,percentChange(ask[0],ask),'r')


    ax1_2 = ax1.twinx()

    ax1_2.fill_between(date, 0, (ask-bid), facecolor='g',alpha=.3)

    plt.subplots_adjust(bottom=.23)

    plt.show()
    