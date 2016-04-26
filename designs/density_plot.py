from __future__ import print_function

try:
    from PyDesign import *
    import PyDesign
except ImportError, e:
    print(e.message)
    exit(1)

import RNAdesign as rd
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

def plot_sequence_objective(args):

    x_obj = []
    y_obj = []
    path = '*.csv'
    if args.directory is not None:
        if not os.path.exists(args.directory):
            print("Path does not exist!", file=sys.stderr)
            exit()
        path = os.path.join(args.directory, path)

    filenames = {}
    if args.file is not None:
        if ";" in args.file:
            filenames = args.file.split(";")
        else:
            filenames.append(args.file)
    else:
        filenames = glob.glob(path)

    where_is_seq_from = {}
    total = 0
    
    for filename in filenames:
        with open(filename) as sampled_seq:
            reader = csv.reader(sampled_seq, delimiter=";")
            header = reader.next()

            # get information of intial sequence
            initial_sequence = reader.next()
            x_center = float(initial_sequence[0])
            y_center = float(initial_sequence[1])


            x_origin = 0 # = x_center-x_center
            y_origin = 0
            total += 1
 
            x_obj.append(x_origin)
            y_obj.append(y_origin)

            #quad3_perc = 0 # NAME
            #which_seq = []
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
    plt.figure()
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
      
    # additional weighting lines 
    if args.more_lines is True:
        #for i in pl.frange(-40, 40, 2):
        for i in pl.frange(-60, 60, 2):
            plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight - i, -ax_lim * args.weight - i], '0.8')
        #plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight - 10, -ax_lim * args.weight - 10], '0.8')
        
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
        '''if args.more_lines is True:  
            step_size = 1    
            limit = 6 -1 #x-achse   
            counts = [0]* (limit+1)'''
        
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
        
    plt.text(ax_lim - 4, ax_lim + 2, quad1, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), ax_lim + 2, quad2, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), -(ax_lim - 2), quad3, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(ax_lim - 4, -(ax_lim - 2), quad4, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})
    
    plt.text(-(ax_lim - 2), ax_lim * args.weight + 0.65, sum_weighted, style='normal', # "hauptgewichtslinie"
    bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
    
    for i in pl.frange(0,7,1):
        plt.text(-10 - i, ax_lim - 2.5, i, style='italic', color='0.4', fontsize=12) # beschriftung "score lines"
    
    for i in range(-20,20,5):   
        plt.text(i, -0.4, '|', style='normal', color='0.4')
  
    plt.pcolormesh(xedges,yedges,Hmasked, vmax = 3000) # color legend
    plt.xlabel('objective 1')
    plt.ylabel('objective 2')
    plt.axhline(color='k')
    plt.axvline(color='k')
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Counts')
    plt.savefig(args.out_file + '.png')
    plt.close()
    
    print('# Output file: %s'% args.out_file + '.png')
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
    if total_datapoints != total:
        print('Number of points in plot and number of sequences are not equal.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot...')
    parser.add_argument("-w", "--weight", type=float, default=0.5**-1, help='Define weighting-factor')
    parser.add_argument("-o", "--out_file", type=str, default= datetime.datetime.now().isoformat(), help='Name graphic')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-a", "--axis_limit", type=int, default=None, help='Axis limit')
    parser.add_argument("-g", "--grid_size", type=int, default=None, help='Grid size')
    parser.add_argument("-c", "--circle_grid", type=float, default=None, help='Circle grid size')
    parser.add_argument("-m", "--more_lines", default=False, action='store_true', help='Additional weighting lines')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.csv format. Separate multiple filenames with ";". If not set reads all *.csv in current directory')
    parser.add_argument("-d", "--directory", type = str, default=None, help='Set directory')
    args = parser.parse_args()

    plot_sequence_objective(args)



