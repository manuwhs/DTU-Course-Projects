import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
import utilities_lib as ul
import matplotlib.gridspec as gridspec

#####  BUILDING FUNCTIONS #####

def savefig(self,file_dir = "./image.png", 
            bbox_inches = 'tight',
            sizeInches = [],  # The size in inches as a list
            close = False,   # If we close the figure once saved
            dpi = 100):      # Density of pixels !! Same image but more cuality ! Pixels
    ## Function to save the current figure in the desired format
    ## Important !! Both dpi and sizeInches affect the number of pixels of the image
    # dpi: It is just for quality, same graph but more quality.
    # sizeInches: You change the proportions of the window. The fontsize and 
    # thickness of lines is the same. So as the graph grows, they look smaller.
    F = self.fig  # ?? NO work ? 
#    F = pylab.gcf()

    Winches, Hinches = F.get_size_inches()  # 8,6
    
#    print Winches, Hinches 
    
    if (len(sizeInches) > 0):
        F.set_size_inches( (sizeInches[0], sizeInches[1]) )
    
    self.fig.savefig(file_dir,
            bbox_inches='tight',
            dpi = dpi
            )
    # Transform back
    F.set_size_inches((Winches, Hinches) )  

    if (close == True):
        plt.close()
        
#    gl.savefig('foodpi50.png', bbox_inches='tight',  dpi = 50)
#    gl.savefig('foodpi100.png', bbox_inches='tight',  dpi = 100)
#    gl.savefig('foodpi150.png', bbox_inches='tight',  dpi = 150)
#
#    gl.savefig('foosize1.png', bbox_inches='tight',  sizeInches = [3,4])
#    gl.savefig('foosize2.png', bbox_inches='tight',  sizeInches = [6,8])
#    gl.savefig('foosize3.png', bbox_inches='tight',  sizeInches = [8,11])
#    gl.savefig('foosize4.png', bbox_inches='tight',  sizeInches = [10,14])

def set_subplots(self, nr, nc, dim = "2D"):
    # State a subplot partitition of a new figure.
    # nr is the number of rows of the partition
    # nc is the numbel of columns
    
    self.subplotting = 1  
    # Variable that indicates that we are subplotting
    # So when nf = 1, if this variable is 0, then 
    # we do not create a new figure but we paint in a
    # different window 
    
    self.nc = nc
    self.nr = nr
    
    self.init_figure(dim = dim)
    
    self.G = gridspec.GridSpec(nr, nc)
    
    self.nci = 0
    self.nri = 0
    
    # We can select the subplot with  plt.subplot(G[r, c])

    self.first_subplot = 1
    
def next_subplot(self, dim = "2D"):
    # Moves to the next subpot to plot.
    # We move from left to right and up down.


    if (self.first_subplot == 1): # If it is the first plot
        # then we do not set the new subplot due to nf = 1
        # This way all the subplottings are the same
        # Select first plot.
        self.first_subplot = 0
    else:
        self.nci = (self.nci + 1) % self.nc
        if (self.nci == 0): # If we started new row
            self.nri = (self.nri + 1) % self.nr
            
        if (self.nri == (self.nr-1) and self.nci == (self.nc-1)): # If we are in the last graph 
            self.subplotting = 0
    # Select next plot.
    if (dim == "2D"):
        self.axes = plt.subplot(self.G[self.nri, self.nci])
    elif(dim == "3D"):
        self.axes = plt.subplot(self.G[self.nri, self.nci], projection='3d')

    if (self.nci + self.nri == 0):
        plt.tight_layout()  # So that the layout is tight

    
def init_figure(self, dim = "2D"):
    # This function initializes the data structures of the new figure
    fig = plt.figure()  
    self.fig = fig
    
    # Get the axes of the plot. 
    # We might want to plot different axes sometimes.
    if (dim == "2D"):
        self.axes = plt.subplot2grid((40,40), (0,0), rowspan=40, colspan=40)  
    elif (dim == "3D"): # No really need, since the 3D func create axis anyway
        self.axes = plt.axes(projection='3d')
    
def set_labels(self, labels, fontsize, loc = 1):
    # This function sets the labels of the graph when created
    # labels: If new figure, we expect 3 strings.
    # Set the main labels !!

    if (len(labels) > 0):
        title = labels[0]
        self.axes.title.set_text(title)
        
    if (len(labels) > 1):
        xlabel = labels[1]
        plt.xlabel(xlabel, fontsize=fontsize)
    if (len(labels) > 2):
        ylabel = labels[2]
        plt.ylabel(ylabel, fontsize=fontsize)

    self.axes.title.set_fontsize(fontsize=fontsize)

def update_legend(self, legend, NcY, loc = 1):
       # If labels specified
    self.axes.legend(loc=loc)
    if(len(legend) > 0):
        self.legend.extend(legend)
    else:
        self.legend.extend(["Line"]*NcY)
    # Plot the legend
#    self.axes.legend(self.legend, loc=loc)
    self.axes.legend(loc=loc)

#    l = plt.legend()
    if (self.axes.legend()):
        self.axes.legend().set_zorder(1000) # Set legend on top
def format_axes (self,nf = 0, fontsize = -1):
    if (len(self.ticklabels) > 0):
#        self.axes.set_xticklabels(self.ticklabels)
        plt.xticks(self.X, self.ticklabels)
        
        self.ticklabels = []  # Delete them
        
    ## Also, rotate the labels
    if (nf == 1):
        for label in self.axes.xaxis.get_ticklabels():
            label.set_rotation(45)
        plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    plt.subplots_adjust(bottom=.15)  # Done to view the dates nicely

    ax = self.axes
    if (fontsize != -1):
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)   
        
def format_plot(self):
    plt.grid()
    plt.show()

def preprocess_data(self,X,Y):
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
    
    # If the type are labels 
    
    if (X.size > 0):
        if (type(X[0,0]).__name__ == "str" or type(X[0,0]).__name__ == "string_"):
            self.ticklabels = X.T.tolist()[0]
            X = ul.fnp([])

    # If we have been given an empty X
    if (X.size == 0):
        # If we are also given an empty Y
        if (Y.size == 0): 
            return -1 # We get out.
            # This could be used to create an empty figure if
            # fig = 1
        # If we wanted to plot something but we did not specify X.
        # We use the previous X, if it does not exist, we use range()
        else:  
            X = ul.fnp(range(NpY)) # Create X axis
                
    # The given X becomes the new axis
    self.X = X    
    self.Y = Y

def preprocess_dates(X):
    # Dealing with dates !
    if (type(X[0,0]).__name__ == "datetime64"):
        X = pd.to_datetime(X).T.tolist()  #  DatetimeIndex
    else:  #  DatetimeIndex
        X = X.T.tolist()[0]  
    
    return X
    
def get_color(self, color):
    # This function outputs the final color to print for a given
    # plotting

    if (type(color).__name__ == "NoneType"):
        # If no color specified. We use one of the list
        colorFinal = self.colors[self.colorIndex]
        self.colorIndex = (self.colorIndex + 1) %len(self.colors)
        return colorFinal
        
def filler(self, X, Yi, color, alpha):
    # This function fills a unique plot.
## We have to fucking deal with dates !!
# The fill function does not work properly with datetime64
        X =  preprocess_dates(X)
        plt.fill_between(X, ul.fnp(Yi).T.tolist()[0], 
                         color = color, alpha = 0.3) 

def figure_management(self, nf, na, labels, fontsize, dim = "2D"):
    # This function si suposed to deal with everything that
    # has to do with initializating figure, axis, subplots...

    if (nf == 1):     
        if (self.subplotting == 1):
            self.next_subplot(dim = dim) # We plot in a new subplot !
        else:
            self.init_figure()  # We create a new figure !!
        self.set_labels(labels, fontsize)
        self.legend = []
        
    # If we want to create a new axis
    elif (na == 1):
        # Define axis to plot volume
        self.axes = self.axes.twinx()  # Create a twin axis

def add_text(self, positionXY = [], text = r'an equation: $E=mc^2$',fontsize = 15):
    
    ## PositonXY should be given in termns of the X and Y axis variables
    if (len(positionXY) == 0):
        positionXY = [0,0]
        
    self.axes.text(positionXY[0], positionXY[1], text, fontsize=fontsize)