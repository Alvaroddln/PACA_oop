import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import ipywidgets as widgets

style = {'description_width': 'initial'}

colors = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange', 'tab:cyan', 'tab:purple', 'tab:pink', 'tab:brown',
          'tab:olive', 'lightcoral', 'seagreen', 'mediumpurple', 'gold', 'sandybrown', 'deepskyblue', 'palegreen']

wells_in_plate = ['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9', 'A10', 'A11', 'A12',
                  'B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12',
                  'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11', 'C12',
                  'D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'D12',
                  'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8', 'E9', 'E10', 'E11', 'E12',
                  'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12',
                  'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9', 'G10', 'G11', 'G12',
                  'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12']

def plot_plate(DF_def, wells_def, time, samples_def, blanks_def):
    max_count = 583 * len(DF_def) #parameter for progress bar
    #progress bar widgets parameters
    f = widgets.IntProgress(min=0, max=max_count,
                    description='Plotting plate datasets... ',
                    bar_style='success',
                    style=style) # instantiate the bar
    display(f) # display the bar
    #cross function
    x = np.linspace(0,1)
    y = x
    y2 = -x+1
    #create plots
    for dataframe in DF_def:
        u = max(DF_def[dataframe].max()) #get y max
        f.value += 1
        d = min(DF_def[dataframe].min()) #get y min
        f.value += 1
        fig, axs = plt.subplots(8, 12, figsize=(15,10)) #create 96 plots in 8x12 array
        f.value += 1
        fig.suptitle(dataframe)
        f.value += 1
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)
        f.value += 1
        #create each plot
        for i, ax in zip(wells_in_plate, axs.ravel()):
            f.value += 1
            #if well is in analysis 
            if i in wells_def:
                ax.plot(time, DF_def[dataframe][i], label=i)
                ax.set_ylim(d, u)
                ax.set_xlim(0,time[-1])
                for sample, c in zip(samples_def, colors):
                  if i in samples_def[sample]:
                    ax.patch.set_facecolor(c)
                    ax.patch.set_alpha(0.2)
                #display blank as text
                if i in blanks_def['blank']:
                    ax.text(time[-1]/2, u/2, 'BLANK', ha='center')
            #else well is outside analysis plot red cross
            else:
                ax.plot(x, y, c='tab:red')
                ax.plot(x, y2, c='tab:red')
                ax.set_ylim(0,1)
                ax.set_xlim(0,1)
                
            f.value += 1
            ax.set_title(i)
            f.value += 1
            ax.set_xticks([])
            f.value += 1
            ax.set_yticks([])
            f.value += 2
        plt.show()
        f.value += 1
        print('')
        f.value += 1# signal to increment the progress bar