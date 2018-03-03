import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
import utilities_lib as ul
import matplotlib.gridspec as gridspec

def Heiken_Ashi_graph(date, data,volume,labels,new_fig = 1):
    r_close = data["Close"].values
    r_open = data["Open"].values
    r_max = data["High"].values
    r_min = data["Low"].values
    
    x_close  = (r_close + r_open + r_max + r_min)/4
    x_open = (r_close[1:] + x_close[:-1])/2  # Mean of the last 
    
    # Insert the first open sin we cannot calcualte the prevoius
    # The other opion is to eliminate the first of everyone
    x_open = np.insert(x_open, 0, r_open[0], axis = 0)
    
#    print x_close.shape, x_open.shape
    
    x_max = np.max(np.array([r_max,x_close,x_open]), axis = 0)
    x_min = np.min(np.array([r_min,x_close,x_open]), axis = 0)
    
    
    
    data = np.array([x_close,x_open,x_max,x_min])
    
    Velero_graph(date, data,volume, labels,new_fig)
    
#    x_close  = np.mean(data,0)
#    x_open = data[0][1:] + upper_box[:-1]  # Mean of the last 
def Velero_graph(date, data, volume, labels,new_fig = 1):
    """ This function plots the Heiken Ashi of the data 
    data[4][Ns]
    -> data[0] = Close
    -> data[1] = Open
    -> data[2] = Max
    -> data[3] = Min
    """
    
    colorFill = "green"  # Gold
    colorBg = "#7fffd4" # Aquamarine
    colorInc = '#FFD700'
    colorDec = "black" 
    
    data = np.array(data)
    data_shape = data.shape
    
    if (date == []):
        date = range(data_shape[1])
        
    Ns = data_shape[1] # Number of samples
    date_indx = np.array(range(data_shape[1]))
    # The width of the bowes is the same for all and the x position is given by
    # the position 

    ## Obtain the 4 parameters of the square
    ## Box parameters

    r_close = data[0]
    r_open = data[1]
    r_max = data[2]
    r_min = data[3]
    
#    x_close  = np.mean(data,0)
#    x_open = data[0][1:] + upper_box[:-1]  # Mean of the last 
    
    """ WE are gonna plot the Velas in one axis and the volume in others """

    ####### Plot all the cubes !!!!
    fig, ax = plt.subplots()
    
    fig.facecolor = colorBg
    
    for i in range(Ns):
        # Calculate upper and lowe value of the box and the sign
        diff = r_close[i] - r_open[i]
#        print diff
        if (diff >= 0):
            low_box = r_open[i]
            sign = colorInc
        else:
            low_box = r_close[i]
            sign = colorDec
        
        # Create the box
        ax.broken_barh([(date_indx[i] + 0.05, 1 - 0.1)],
                        (low_box, abs(diff)),
                         facecolors=sign)
                         
        # Create the box upper line                
        ax.broken_barh([(date_indx[i] + 0.45, 0.1)],
                        (low_box + abs(diff), r_max[i] - low_box + abs(diff)),
                         facecolors= "red")
                         
        # Create the box lower line                
        ax.broken_barh([(date_indx[i] + 0.45, 0.1)],
                        (r_min[i] , low_box - r_min[i]),
                         facecolors= "red")
    
    """ PLOT VOLUME """ 
    ax1_2 = ax.twinx()
    
    #ax1_2.plot(date, (ask-bid))
#    print volume.shape
#    print len(date)
    ax1_2.bar(date, volume, facecolor= colorFill,alpha=.5)
#    ax1_2.fill_between(date, 0, volume, facecolor= colorFill,alpha=.5)
    
#    broken_barh ( xranges, yrange, **kwargs)
#        xranges	sequence of (xmin, xwidth)
#        yrange	 sequence of (ymin, ywidth)
        
    plt.title(labels[0])
    plt.xlabel(labels[1])
    plt.ylabel(labels[2])
    
    if (len(labels) >3 ):
        plt.legend(labels[3])
    
#    plt.grid()
    plt.show()
    

def TrendVelero_graph(date, data, volume, labels,new_fig = 1):
    """ This function plots the Heiken Ashi of the data 
    data[4][Ns]
    -> data[0] = Max
    -> data[1] = Min
    """
    
    colorFill = "green"  # Gold
    colorBg = "#7fffd4" # Aquamarine
    
    colorInc = '#FFD700'
    colorDec = "black" 
    
    colorDiv = 'red'  # Divergence
    colorCon = "blue"    # Convergence
    
    data = np.array(data)
    data_shape = data.shape
    
    if (date == []):
        date = range(data_shape[1])
        
    Ns = data_shape[1] # Number of samples
    date_indx = np.array(range(data_shape[1]))
    # The width of the bowes is the same for all and the x position is given by
    # the position 

    ## Obtain the 4 parameters of the square
    ## Box parameters

    r_max = data[0]
    r_min = data[1]
    
#    x_close  = np.mean(data,0)
#    x_open = data[0][1:] + upper_box[:-1]  # Mean of the last 
    
    """ WE are gonna plot the Velas in one axis and the volume in others """

    ####### Plot all the cubes !!!!
    fig, ax = plt.subplots()
    
    fig.facecolor = colorBg
    
    for i in range(Ns):
        # Calculate upper and lowe value of the box and the sign
        diff_max = r_max[i] - r_max[i-1]
        diff_min = r_min[i] - r_min[i-1]
        low_box = r_min[i]
        diff = r_max[i] - r_min[i]
        
#        print diff
        if (diff_max >= 0)&(diff_min >= 0):  # Increasing
            sign = colorInc
            
        elif(diff_max < 0)&(diff_min < 0):   # Decreasing
            sign = colorDec
            
        elif(diff_max > 0)&(diff_min < 0):   # More Variance 
            sign = colorDiv
            
        elif(diff_max < 0)&(diff_min > 0):   # Less variance
            sign = colorCon  
            
        # Create the box
        ax.broken_barh([(date_indx[i] + 0.05, 1 - 0.1)],
                        (low_box, abs(diff)),
                         facecolors=sign)
                         
    """ PLOT VOLUME """ 
    ax1_2 = ax.twinx()
    
    #ax1_2.plot(date, (ask-bid))
#    print volume.shape
#    print len(date)
    ax1_2.bar(date, volume, facecolor= colorFill,alpha=.5)
#    ax1_2.fill_between(date, 0, volume, facecolor= colorFill,alpha=.5)
    
#    broken_barh ( xranges, yrange, **kwargs)
#        xranges	sequence of (xmin, xwidth)
#        yrange	 sequence of (ymin, ywidth)
        
    plt.title(labels[0])
    plt.xlabel(labels[1])
    plt.ylabel(labels[2])
    
    if (len(labels) >3 ):
        plt.legend(labels[3])
    
#    plt.grid()
    plt.show()
    

def scatter_graph(x,y,labels,new_fig = 1):
    y = np.array(y)
#    print y.shape
#    print x.shape
    if (new_fig == 1):
        plt.figure(figsize = [w,h])
        
    if (x == []): # If x is empty
        print "X was empty"
        y_shape = y.shape
        x = np.array(range(y_shape[0]))
        
    print x.shape, y.shape
        # TODO for some reason doesnt let me print matrix so convert to list
    
    plt.scatter(x.tolist(),y.tolist(), color = np.random.uniform(0,1,(1,3)), lw = 3, alpha = 0.2)

    plt.title(labels[0])
    plt.xlabel(labels[1])
    plt.ylabel(labels[2])
    
    if (len(labels) > 3 ):
        plt.legend(labels[3])
    
    plt.grid()
    plt.show()
    

