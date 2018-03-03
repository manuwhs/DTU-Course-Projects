import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import utilities_lib as ul

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def preproces_data_3D(self, X,Y):
    
   # Preprocess the variables X and Y
   ### X axis date properties, just to be seen clearly
   ### First we transform everything to python format
    X = ul.fnp(X)
    Y = ul.fnp(Y)

    # Each plot can perform several plotting on the same X axis
    # So there could be more than one NcX and NcY
    # NpX and NpY must be the same.
    NpX, NcX = X.shape
    NpY, NcY = Y.shape
    
    if (X.size > 0):
        if (type(X[0,0]).__name__ == "str" or type(X[0,0]).__name__ == "string_"):
            self.Xticklabels = X.T.tolist()[0]
            X = ul.fnp(range(X.size))
            
    if (Y.size > 0):
        if (type(Y[0,0]).__name__ == "str" or type(Y[0,0]).__name__ == "string_"):
            self.Yticklabels = Y.T.tolist()[0]
            Y = ul.fnp(range(X.size)) 
            
    # The given X becomes the new axis
    self.X = X    
    self.Y = Y
            
    return X,Y

def format_axes_3D (self,nf = 0, fontsize = -1):
    
    if (len(self.Xticklabels) > 0):
#        self.axes.set_xticklabels(self.ticklabels)
        # + 0.5 to center the tags
        plt.xticks(self.X + 0.5, self.Xticklabels)
        
        self.Xticklabels = []  # Delete them
    
    if (len(self.Yticklabels) > 0):
#        self.axes.set_xticklabels(self.ticklabels)
    
        plt.yticks(self.Y + 0.5, self.Yticklabels)
        
        self.Yticklabels = []  # Delete them

    ax1 = self.axes
    ax = ax1
    [t.set_va('center') for t in ax1.get_yticklabels()]
    [t.set_ha('left') for t in ax1.get_yticklabels()]
    [t.set_va('center') for t in ax1.get_xticklabels()]
    [t.set_ha('right') for t in ax1.get_xticklabels()]
    [t.set_va('center') for t in ax1.get_zticklabels()]
    [t.set_ha('left') for t in ax1.get_zticklabels()]
    
    ## Also, rotate the labels
#    if (nf == 1):
#        for label in self.axes.xaxis.get_ticklabels():
#            label.set_rotation(45)
#        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
#    plt.subplots_adjust(bottom=.15)  # Done to view the dates nicely

    if (fontsize != -1):
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)   
        for tick in ax.zaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)     
            
def plot_3D (self, Xgrid,Ygrid, Zvalues,
        labels = [],     # Labels for tittle, axis and so on.
        legend = [],     # Legend for the curves plotted
        nf = 1,          # New figure
        na = 0,          # New axis. To plot in a new axis
        # Basic parameters that we can usually find in a plot
        color = None,    # Color
        lw = 2,          # Line width
        alpha = 1.0,      # Alpha
        
        fontsize = 15,   # The font for the labels in the title
        fontsizeL = 10,  # The font for the labels in the legeng
        fontsizeA = 20,  # The font for the labels in the axis
        loc = 1,
        ):
        # Plots the training and Validation score of a realization

    X = np.array(Xgrid)
    Y = np.array(Ygrid)
    X, Y = np.meshgrid(X, Y)
    
    Z = ul.convert_to_matrix(Zvalues).T

    # Management of the figure
    self.figure_management(nf, na, labels, fontsize, dim = "3D")
    
#    self.ax = self.fig.gca(projection='3d')

    fig = self.fig
    ax = self.axes
    
    if (nf == 1):
        color = 'copper'
    else:
        color = cm.coolwarm
        
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, 
                           cmap=color,      #  cm.coolwarm   ('Sequential',  ['Greys'])
                
                           linewidth=0, antialiased=False)
    
    ax.set_zlim(np.min(Z.flatten()), np.max(Z.flatten()))
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    ax.set_zlabel('Z Label')
    
#    fig.colorbar(surf, shrink=0.5, aspect=5)
 
    plt.show()
    
    
def bar_3D (self, Xgrid,Ygrid, Zvalues,
        labels = [],     # Labels for tittle, axis and so on.
        legend = [],     # Legend for the curves plotted
        nf = 1,          # New figure
        na = 0,          # New axis. To plot in a new axis
        # Basic parameters that we can usually find in a plot
        color = None,    # Color
        lw = 2,          # Line width
        alpha = 1.0,      # Alpha
        
        width = 1.0,
        fontsize = 15,   # The font for the labels in the title
        fontsizeL = 10,  # The font for the labels in the legeng
        fontsizeA = 20,  # The font for the labels in the axis
        loc = 1,
        ):
        # Plots the training and Validation score of a realization

    # Initial position of the cube in 3D
    X = np.array(Xgrid)
    Y = np.array(Ygrid)
    
    self.preproces_data_3D(X,Y)
    
    X = self.X
    Y = self.Y
    
    xpos, ypos = np.meshgrid(X, Y)
    
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(ul.fnp(Zvalues).shape).flatten()
    
#    print xpos.shape, ypos.shape, zpos.shape

    ## Lenght of the cubes
    dx = np.ones(xpos.size).flatten() * width
    dy = np.ones(ypos.size).flatten() * width
    dz = ul.fnp(Zvalues).flatten()

    self.figure_management(nf, na, labels, fontsize, dim = "3D")
    
    fig = self.fig
    ax = self.axes
    
    if (nf == 1):
        color = 'copper'
    else:
        color = cm.coolwarm
        
    ax.bar3d(xpos, ypos, zpos, 
             dx, dy, dz, 
             cmap="copper",  color='#00ceaa', alpha = alpha) 
 
    plt.show()
    self.format_axes_3D(nf, fontsize = fontsizeA)
    


