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

            x_obj.append(x_origin)
            y_obj.append(y_origin)

            #quad3_perc = 0 # NAME
            #which_seq = []
            for line in reader:
                x_new = float(line[0])
                y_new = float(line[1])

                # calculate relative positions to initial sequence
                x = x_new - x_center
                y = y_new - y_center
                
                x_obj.append(x)
                y_obj.append(y)
                
                # information for entries in quadrant 3
                '''if x < 0 and y < 0: # in quadrant 3
                    score = float(line[2])
                    seq = line[3] 
                    quad3_perc += 1 # only quadrant 3 
                    which_seq.append(seq)
                    which_seq.append(score) 
                    which_seq.append(x)
                    which_seq.append(x_new)
                    which_seq.append(y)
                    which_seq.append(y_new)
                            
            if quad3_perc != 0:   
                which_seq.insert(0, quad3_perc)    
                where_is_seq_from[filename] = which_seq'''

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
  
    # plot "weighting lines"    
    plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight, -ax_lim * args.weight], '0.8') # through origin
    plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight - 2.5, -ax_lim * args.weight - 2.5], '0.8') 
    
    # additional weighting lines for 5000 mutations
    if args.more_lines is True:
        for i in pl.frange(0.5, 2.5, 0.5):
            plt.plot([-ax_lim, ax_lim], [ax_lim * args.weight - i, -ax_lim * args.weight - i], '0.8')
        
       
    if args.grid_size is not None:
        grid = args.grid_size
        for i in range(grid, ax_lim, grid): 
            plt.plot([-ax_lim, ax_lim], [i, i], '0.75')
            plt.plot([-ax_lim, ax_lim], [-i, -i], '0.75')
            
        for i in range(grid, ax_lim,grid): 
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
    
    # counts per "additional lines"
    up_quarter_line = 0
    d_quarter_line = 0
    
    # counts for weighting lines
    weighted1 = 0
    weighted2 = 0
    
    quart0 = 0
    quart1 = 0
    quart2 = 0
    quart3 = 0
    quart4 = 0
        
    for i in range(0, len(x_obj)):
        if x_obj[i] > 0 and y_obj[i] > 0:
            quad1 += 1
            
        if x_obj[i] < 0 and y_obj[i] > 0:
            quad2 += 1
            if y_obj[i] < (-args.weight * x_obj[i]):
                weighted1 += 1
                if y_obj[i] >= (-args.weight * x_obj[i] - 2.5): 
                    up_quarter_line += 1
                else:
                    d_quarter_line += 1
                
        if x_obj[i] < 0 and y_obj[i] < 0:
            quad3 += 1
            if y_obj[i] >= (-args.weight * x_obj[i] - 2.5): 
                up_quarter_line += 1
            else:
                d_quarter_line += 1
                
        if x_obj[i] > 0 and y_obj[i] < 0:
            quad4 += 1
            if y_obj[i] < (-args.weight * x_obj[i]):
                weighted2 += 1
                if y_obj[i] >= (-args.weight * x_obj[i] - 2.5): 
                    up_quarter_line += 1
                else:
                    d_quarter_line += 1

    if args.more_lines is True:
        for i in range(0, len(x_obj)):           
            if x_obj[i] < 0 and y_obj[i] > 0:
                #quad2 += 1
                if y_obj[i] < (-args.weight * x_obj[i]): # unter weightingline
                    #weighted1 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 2.5) and y_obj[i] < (-args.weight * x_obj[i] - 2): 
                        quart4 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 2) and y_obj[i] < (-args.weight * x_obj[i] - 1.5):
                        quart3 +=1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 1.5) and y_obj[i] < (-args.weight * x_obj[i] - 1):
                        quart2 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 1) and y_obj[i] < (-args.weight * x_obj[i] - 0.5):
                        quart1 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 0.5) and y_obj[i] < (-args.weight * x_obj[i]):
                        quart0 += 1
                   
            if x_obj[i] < 0 and y_obj[i] < 0:
                #quad3 += 1
                if y_obj[i] >= (-args.weight * x_obj[i] - 2.5) and y_obj[i] < (-args.weight * x_obj[i] - 2): 
                    quart4 += 1
                if y_obj[i] >= (-args.weight * x_obj[i] - 2) and y_obj[i] < (-args.weight * x_obj[i] - 1.5):
                    quart3 +=1
                if y_obj[i] >= (-args.weight * x_obj[i] - 1.5) and y_obj[i] < (-args.weight * x_obj[i] - 1):
                    quart2 += 1
                if y_obj[i] >= (-args.weight * x_obj[i] - 1) and y_obj[i] < (-args.weight * x_obj[i] - 0.5):
                    quart1 += 1
                if y_obj[i] >= (-args.weight * x_obj[i] - 0.5) and y_obj[i] < (-args.weight * x_obj[i]):
                    quart0 += 1
                
                    
            if x_obj[i] > 0 and y_obj[i] < 0:
                #quad4 += 1
                if y_obj[i] < (-args.weight * x_obj[i]):
                    #weighted2 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 2.5) and y_obj[i] < (-args.weight * x_obj[i] - 2): 
                        quart4 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 2) and y_obj[i] < (-args.weight * x_obj[i] - 1.5):
                        quart3 +=1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 1.5) and y_obj[i] < (-args.weight * x_obj[i] - 1):
                        quart2 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 1) and y_obj[i] < (-args.weight * x_obj[i] - 0.5):
                        quart1 += 1
                    if y_obj[i] >= (-args.weight * x_obj[i] - 0.5) and y_obj[i] < (-args.weight * x_obj[i]):
                        quart0 += 1


    sum_weighted = weighted1 + weighted2 + quad3
    
    for entry in where_is_seq_from:
        x = where_is_seq_from[entry]
        num_in_3 = float(x[0])
        percentage = num_in_3/quad3
    #    print(entry) # entry is filename
        x.insert(0, percentage) #replace count with percentage
        where_is_seq_from[entry] = x
    #   print(where_is_seq_from[entry])
        
  
        

    plt.text(ax_lim - 4, ax_lim + 2, quad1, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), ax_lim + 2, quad2, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), -(ax_lim - 2), quad3, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(ax_lim - 4, -(ax_lim - 2), quad4, style='italic',
    bbox={'facecolor':'white', 'alpha':0.5, 'pad':10})

    plt.text(-(ax_lim - 2), ax_lim * args.weight + 0.75, sum_weighted, style='italic', # "hauptgewichtslinie"
    bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
    
    plt.text(-(ax_lim - 2), ax_lim * args.weight - 2.5, up_quarter_line, style='italic', # zwischen 2.5 und hauptgewicht
    bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})
    
    plt.text(-(ax_lim - 2), -3, d_quarter_line, style='italic', # "unter zweiter gewichtslinie"
    bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

    plt.pcolormesh(xedges,yedges,Hmasked)
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
        print('# 0: %s'% quart0)
        print('# 1: %s'% quart1)
        print('# 2: %s'% quart2)
        print('# 3: %s'% quart3)
        print('# 4: %s'% quart4)
    
    #print(where_is_seq_from)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot...')
    parser.add_argument("-w", "--weight", type=float, default=0.5, help='Define weighting-factor')
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



