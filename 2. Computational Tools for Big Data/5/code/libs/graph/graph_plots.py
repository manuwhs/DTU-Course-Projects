import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
import utilities_lib as ul
import matplotlib.gridspec as gridspec

import graph_basic as grba

def plot(self, X = [],Y = [],  # X-Y points in the graph.
        labels = [],     # Labels for tittle, axis and so on.
        legend = [],     # Legend for the curves plotted
        nf = 1,          # New figure
        na = 0,          # New axis. To plot in a new axis
        # Basic parameters that we can usually find in a plot
        color = None,    # Color
        lw = 2,          # Line width
        alpha = 1.0,      # Alpha
        
        fontsize = 20,   # The font for the labels in the title
        fontsizeL = 10,  # The font for the labels in the legeng
        fontsizeA = 15,  # The font for the labels in the axis
        loc = 1,
        
        ### Super Special Shit !!
        fill = 0  #  0 = No fill, 1 = Fill and line, 2 = Only fill
        ):         

    ## Preprocess the data given so that it meets the right format
    self.preprocess_data(X,Y)
    X = self.X
    Y = self.Y
    
    NpX, NcX = X.shape
    NpY, NcY = Y.shape
    
    # Management of the figure
    self.figure_management(nf, na, labels, fontsize)

    ##################################################################
    ############### CALL PLOTTING FUNCTION ###########################
    ##################################################################
    ## TODO. Second case where NcY = NcX !!

    for i in range(NcY):  # We plot once for every line to plot
        self.zorder = self.zorder + 1
        colorFinal = self.get_color(color)
        if (i >= len(legend)):
            plt.plot(X,Y[:,i], lw = lw, alpha = alpha, 
                     color = colorFinal, zorder = self.zorder)
        else:
#            print X.shape
#            print Y[:,i].shape
            plt.plot(X,Y[:,i], lw = lw, alpha = alpha, color = colorFinal,
                     label = legend[i], zorder = self.zorder)
        
        if (fill == 1):  ## Fill this shit !!
            self.filler(X,Y[:,i],colorFinal,alpha)
    
    self.update_legend(legend,NcY,loc = loc)
    self.format_axes(nf, fontsize = fontsizeA)
    

    if (na == 1 or nf == 1):
        self.format_plot()
    
    return 0

def scatter(self, X = [],Y = [],  # X-Y points in the graph.
        labels = [],     # Labels for tittle, axis and so on.
        legend = [],     # Legend for the curves plotted
        nf = 1,          # New figure
        na = 0,          # New axis. To plot in a new axis
        # Basic parameters that we can usually find in a plot
        color = None,    # Color
        lw = 2,          # Line width
        alpha = 1.0,      # Alpha
        
        fontsize = 20,   # The font for the labels in the title
        fontsizeL = 10,  # The font for the labels in the legeng
        fontsizeA = 15,  # The font for the labels in the axis
        loc = 1
        ):         

    ## Preprocess the data given so that it meets the right format
    self.preprocess_data(X,Y)
    X = self.X
    Y = self.Y
    
    NpX, NcX = X.shape
    NpY, NcY = Y.shape
    
    # Management of the figure
    self.figure_management(nf, na, labels, fontsize)

    ##################################################################
    ############### CALL SCATTERING FUNCTION ###########################
    ##################################################################
    
    ## We asume that X and Y have the same dimensions
    colorFinal = self.get_color(color)
    X =  grba.preprocess_dates(X)
    
    self.zorder = self.zorder + 1
    if (len(legend) == 0):
        plt.scatter(X,Y, lw = lw, alpha = alpha, 
                    color = colorFinal, zorder = self.zorder)
    else:
        plt.scatter(X,Y, lw = lw, alpha = alpha, color = colorFinal,
                    label = legend[0], zorder = self.zorder)
    self.format_axes(nf, fontsize = fontsizeA)
    self.update_legend(legend,NcY,loc = loc)

    if (na == 1 or nf == 1):
        self.format_plot()
    
    return 0

def get_barwidth(X):
    # The Xaxis could be dates and so on, so we want to calculate
    # the with of this bastard independently of that

    if (type(X[0]).__name__ == "Timestamp"):
        width = (X[1] - X[0]) 
        width = width.days 
    else:
    # Set the width so that it will fill all the space
        width = (X[1] - X[0]) #/((10**9)*60*60*24)
#        print width
        
    width = float(width)
    
    return width
    
def bar(self, X = [],Y = [],  # X-Y points in the graph.
        labels = [],     # Labels for tittle, axis and so on.
        legend = [],     # Legend for the curves plotted
        nf = 1,          # New figure
        na = 0,          # New axis. To plot in a new axis
        # Basic parameters that we can usually find in a plot
        color = None,    # Color
        width = -1.0,       # Rectangle width
        alpha = 1.0,     # Alpha
        
        fontsize = 20,   # The font for the labels in the title
        fontsizeL = 10,  # The font for the labels in the legeng
        fontsizeA = 15,  # The font for the labels in the axis
        loc = 1
        ):         

    ## Preprocess the data given so that it meets the right format
    self.preprocess_data(X,Y)
    X = self.X
    Y = self.Y
    
    NpX, NcX = X.shape
    NpY, NcY = Y.shape
    
#    print X.shape, Y.shape
    # Management of the figure
    self.figure_management(nf, na, labels, fontsize)

    ##################################################################
    ############### CALL SCATTERING FUNCTION ###########################
    ##################################################################
    
    ## We asume that X and Y have the same dimensions
    colorFinal = self.get_color(color)
    
    X =  grba.preprocess_dates(X)
    
    if (width < 0):
        width = get_barwidth(X)
    
    if (len(legend) == 0):
        self.axes.bar(X, Y, width = width, align='center',
                      facecolor= colorFinal,alpha=alpha)
    else:
        self.axes.bar(X, Y, width = width, align='center',
                      facecolor= colorFinal,alpha=alpha,
                      label = legend[0])
    self.format_axes(nf, fontsize = fontsizeA)
    self.update_legend(legend,NcY,loc = loc)

    # When nf = 0 and na = 0, we lose the grid for some reason.
    if (na == 1 or nf == 1):
        self.format_plot()
    
    return 0
    
    
def plot2(self):
    fig, ax1 = plt.subplots(figsize = [self.w,self.h])   # Get the axis and we duplicate it as desired
    self.figure = fig
    axis = [ax1]
    print self.labels
    for i in range(self.nplots):
        if (i > 0):
            axis.append(ax1.twinx())  # Create a duplicate of the axises (like new superpeusto plot)
    
        axis[i].plot(self.plot_x[i],self.plot_y[i],lw = self.lw)
        axis[i].set_xlabel('time (s)')
        axis[i].set_ylabel('exp')
        axis[i].legend (self.labels[i])
#            for tl in axis[i].get_yticklabels():
#                tl.set_color('b')
    plt.grid()
    plt.show()

