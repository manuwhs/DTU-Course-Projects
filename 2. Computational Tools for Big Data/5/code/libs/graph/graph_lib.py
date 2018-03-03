import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import utilities_lib as ul

import graph_basic as grba
import graph_plots as grpl
import graph_advanced as grad
import graph_3D as gr3D
w = 20  # Width of the images
h = 12   # Height of the images


class CGraph ():
    
    def __init__(self,w = 20, h = 12, lw = 2):
        self.w = w;    # X-width
        self.h = h;    # Y-height
        self.lw = lw   # line width
        
        self.nplots = 0;  # Number of plots
        self.labels = []  # Labels for the plots
        self.plot_y = []  # Vectors of values to plot x-axis
        self.plot_x = []  # Vectors of values to plot y-axis
        
        # Set of nice colors to iterate over when we do not specify color
        self.colors = ["b", "g","k", "r", "c", "m","y"] 
        self.colorIndex = 0;  # Index of the color we are using.
        self.X = ul.fnp([])
        self.legend = []
        
        self.subplotting = 0;  # In the beggining we are not subplotting
        self.ticklabels = []
        self.Xticklabels = []  # For 3D
        self.Yticklabels = []
        self.zorder = 1   # Zorder for plotting
    
    savefig = grba.savefig
    set_labels = grba.set_labels
    init_figure = grba.init_figure
    update_legend = grba.update_legend
    format_axes = grba.format_axes
    format_plot = grba.format_plot
    preprocess_data = grba.preprocess_data
    get_color = grba.get_color
    figure_management = grba.figure_management
    filler = grba.filler
    
    set_subplots = grba.set_subplots
    next_subplot = grba.next_subplot
    
    plot = grpl.plot
    scatter = grpl.scatter
    bar = grpl.bar
    
    preproces_data_3D = gr3D.preproces_data_3D
    format_axes_3D = gr3D.format_axes_3D
    plot_3D = gr3D.plot_3D
    bar_3D = gr3D.bar_3D
    ###### 

    add_text = grba.add_text
gl = CGraph()
