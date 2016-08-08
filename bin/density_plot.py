from __future__ import print_function

try:
    from PyDesign import *
    import PyDesign
except ImportError, e:
    print(e.message)
    exit(1)

import RNAblueprint as rbp
import argparse
import sys
import time
import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import glob
import os
import pylab as pl
from decimal import Decimal
import math
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

def plot_sequence_objective(args):

    x_obj = []
    y_obj = []
    path = '*.csv'
    if args.directory is not None:
        if not os.path.exists(args.directory):
            print("Path does not exist!", file=sys.stderr)
            exit()
        path = os.path.join(args.directory, path)

    filenames = []
    if args.file is not None:
        if ";" in args.file:
            filenames = args.file.split(";")
        else:
            filenames.append(args.file)
    else:
        filenames = glob.glob(path)

    total = 0
    start_seq = 0
    
    for filename in filenames:
        with open(filename) as sampled_seq:
            reader = csv.reader(sampled_seq, delimiter=";")
            header = reader.next()

            # get information of intial sequence
            initial_sequence = reader.next()
            x_center = float(initial_sequence[0])
            y_center = float(initial_sequence[1])


            #x_origin = 0 # = x_center-x_center
            #y_origin = 0
            start_seq += 1
            total += 1
 
            #x_obj.append(x_origin)
            #y_obj.append(y_origin)

            for line in reader:
            
                x_new = float(line[0])
                y_new = float(line[1])
                
               

                # calculate relative positions to initial sequence
                x = float(x_new - x_center)
                y = float(y_new - y_center)
                total += 1
               
                x_obj.append(x)
                y_obj.append(y)
                
    # find largest x/y and use as limit for both axes
    if max(x_obj) > max(y_obj):
        ax_lim = abs(max(x_obj))
    else:
        ax_lim = abs(max(y_obj))

    # density plot https://oceanpython.org/2013/02/25/2d-histogram/
    # Estimate the 2D histogram
    nbins = 50 # TODO: how to choose?
    
    H, xedges, yedges = np.histogram2d(x_obj, y_obj, bins=nbins)

    H = np.rot90(H)
    H = np.flipud(H)

    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero

    # Plot 2D histogram
    fig = plt.figure()
    
    # grid options
    # lines
    if args.axis_limit is None:
        ax_lim += 1  # +1 to avoid points on axes
        plt.xlim([-ax_lim, ax_lim])
        plt.ylim([-ax_lim, ax_lim])

    else:
        ax_lim = args.axis_limit
        plt.xlim([-ax_lim, ax_lim])
        plt.ylim([-ax_lim, ax_lim])
        
    # additional axes    
    ax1=fig.add_subplot(1,1,1)
    #ax1.tick_params(labeltop=True, labelright=True)
    
    ax1.get_yaxis().set_tick_params(which='both', direction='out')
    ax1.get_xaxis().set_tick_params(which='both', direction='out')
    minor_ticks = np.arange(-30,30,5)
    ax1.set_xticks(minor_ticks, minor = True)
    ax1.set_yticks(minor_ticks, minor = True)
      
    # additional weighting lines 
    if args.more_lines is True:
        #for i in pl.frange(-40, 40, 2):
        for i in pl.frange(-90, 90, 2):
            if i%5 == 0:
                plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight - i, -ax_lim * args.weight - i], '0.7')
            else:
                plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight - i, -ax_lim * args.weight - i], '0.9')
        
    # plot "zero line"   
    plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight, -ax_lim * args.weight], 'm')   
       
    if args.grid_size is not None:
        grid = args.grid_size
        for i in pl.frange(grid, ax_lim, grid): 
            plt.plot([-ax_lim, ax_lim], [i, i], '0.75')
            plt.plot([-ax_lim, ax_lim], [-i, -i], '0.75')
            
        for i in pl.frange(grid, ax_lim,grid): 
            plt.plot( [i, i], [-ax_lim, ax_lim],'0.75')
            plt.plot([-i, -i],[-ax_lim, ax_lim], '0.75')
            
    # circles
    if args.circle_grid is not None:
        ax=fig.add_subplot(1,1,1)
        circle_grid = args.circle_grid
        for i in pl.frange(circle_grid, ax_lim, circle_grid):
            circ =  plt.Circle((0,0),i , color='0.8', fill=False)
            ax.add_patch(circ)

    # counts per quadrants
    quad1 = 0
    quad2 = 0
    quad3 = 0
    quad4 = 0
   
    # counts for weighting lines
    weighted1 = 0
    weighted2 = 0
       
    l = 0 
    if args.more_lines is True:  
            step_size = 1    
            limit = 18 -1 #x-achse   
            counts = [0]* (limit+1)
            
    for i in range(0, len(x_obj)):
        no_add = True
        
        # Quadrant 1
        if x_obj[i] >= 0 and y_obj[i] >= 0:
            quad1 += 1
            no_add = False
            
        # Quadrant 2    
        if x_obj[i] < 0 and y_obj[i] >= 0:
            quad2 += 1
            no_add = False
            if y_obj[i] < (-args.weight * x_obj[i]): # alles unter 0-Linie
                weighted1 += 1
                
                if args.more_lines is True: # zusaetzliche counts
                    for x in pl.frange(0, limit, step_size):
                        if y_obj[i] >= (-args.weight * x_obj[i] - ((x+step_size)*2)) and y_obj[i] < (-args.weight * x_obj[i]-x*2):
                            counts[x] += 1 
                
        # Quadrant 3        
        if x_obj[i] < 0 and y_obj[i] < 0: # alles unter 0-Linie
            quad3 += 1
            no_add = False
            if args.more_lines is True:
                for x in pl.frange(0, limit, step_size):
                        if y_obj[i] >= (-args.weight * x_obj[i] - ((x+step_size)*2)) and y_obj[i] < (-args.weight * x_obj[i]-x*2):
                            counts[x] += 1 
                
        # Quadrant 4      
        if x_obj[i] >= 0 and y_obj[i] < 0:
            quad4 += 1
            no_add = False
            if y_obj[i] < (-args.weight * x_obj[i]): # alles unter 0-Linie
                weighted2 += 1
                if args.more_lines is True:
                    for x in pl.frange(0, limit, step_size):
                        if y_obj[i] >= (-args.weight * x_obj[i] - ((x+step_size)*2)) and y_obj[i] < (-args.weight * x_obj[i]-x*2):
                            counts[x] += 1 
                
                    
        if no_add is True:
            print(str(i) + " " + str(x_obj[i]) +" " + str(y_obj[i]))
            l+=1
            print(l)
   

    sum_weighted = weighted1 + weighted2 + quad3
    
    
    quad1_f = format(quad1, "6,d").replace(",", ".")    
    quad2_f= format(quad2, "6,d").replace(",", ".")
    quad3_f = format(quad3, "6,d").replace(",", ".")
    quad4_f = format(quad4, "6,d").replace(",", ".")
       
    sum_weighted_f = format(sum_weighted, "6,d").replace(",", ".")

    plt.text(-(ax_lim - 2), ax_lim - 9, sum_weighted_f, style='italic',
    bbox={'facecolor':'white', 'edgecolor':'darkmagenta', 'pad':10})
    
    q1 = AnchoredText(quad1_f, prop={'size':14}, frameon=True,loc=1)
    q2 = AnchoredText(quad2_f, prop={'size':14}, frameon=True,loc=2)
    q3 = AnchoredText(quad3_f, prop={'size':14}, frameon=True,loc=3)
    q4 = AnchoredText(quad4_f, prop={'size':14}, frameon=True,loc=4)
    ax1.add_artist(q1)
    ax1.add_artist(q2)
    ax1.add_artist(q3)
    ax1.add_artist(q4)
      
    plt.pcolormesh(xedges,yedges,Hmasked, vmax = 3000) # color legend
    #plt.xlabel('objective 1')
    #plt.ylabel('objective 2')
    plt.axhline(color='k')
    plt.axvline(color='k')
    
    neighbors = total - start_seq 
    neighbors_f = format(neighbors, "6,d").replace(",", ".")
    plt.title('Number of neighbors: %s' % neighbors_f)
    
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.savefig(args.out_file + '.svg')
    plt.close()
    
    print('# Output file: %s'% args.out_file + '.svg')
    print('# Quadrant 1: %s'% quad1)
    print('# Quadrant 2: %s'% quad2)
    print('# Quadrant 3: %s'% quad3)
    print('# Quadrant 4: %s'% quad4)
    
    print('# Quadrant 3++: %s'% sum_weighted)

    if args.more_lines is True:
        print('# Counts per score improvement:')
            
        for i in range(0,len(counts)):
            print('# %s' %i + '-' + str(i+step_size) +': %s'% counts[i])
                   
    total_datapoints = quad1+quad2+quad3+quad4
    print('Total data points in plot: %s' % total_datapoints)
    print('Number of sequences: %s' % total)
    print('# Number of neighbors: %s'% neighbors)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot...')
    parser.add_argument("-w", "--weight", type=float, default=0.5**-1, help='Define weighting-factor')
    parser.add_argument("-o", "--out_file", type=str, default= datetime.datetime.now().isoformat(), help='Name graphic')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-a", "--axis_limit", type=int, default=None, help='Axis limit')
    parser.add_argument("-g", "--grid_size", type=int, default=None, help='Grid size')
    #parser.add_argument("-n", "--neighbors", type=str, default= None, help='Number of neighbors')
    parser.add_argument("-c", "--circle_grid", type=float, default=None, help='Circle grid size')
    parser.add_argument("-m", "--more_lines", default=False, action='store_true', help='Additional weighting lines')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.csv format. Separate multiple filenames with ";". If not set reads all *.csv in current directory')
    parser.add_argument("-d", "--directory", type = str, default=None, help='Set directory')
    args = parser.parse_args()

    plot_sequence_objective(args)



