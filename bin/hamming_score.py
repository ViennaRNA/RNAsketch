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
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import matplotlib.gridspec as gridspec

def box_plot(args):

    
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
    
    data = {}
    distances = []
    
    for filename in filenames:
        with open(filename) as sampled_seq:
            reader = csv.DictReader(sampled_seq, delimiter=";")

            # get information of intial sequence
            initial_sequence = reader.next()
            score_center = float(initial_sequence['score'])
            hamming_center = int(initial_sequence['hamming distance'])
            total += 1

            for line in reader:
                
                score_new = float(line['score'])
                hamming = int(line['hamming distance'])
                
                # calculate relative positions to initial sequence
                rel_score = float(score_new - score_center)
                
                if hamming not in data:
                    data[hamming] = []
                    
                data[hamming].append(rel_score)
                distances.append(hamming)
                total += 1
    
    data_list = []
    
    for keys in data:
        data_list.append(data[keys])
        
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])                                                                
    fig = plt.figure()
            
    #ax1=plt.add_subplot(2,1,2)
    ax1=plt.subplot(gs[1])
    bp = ax1.boxplot(data_list, patch_artist=True)
    plt.axhline(color='0.7')
    ax1.tick_params(labeltop=True)
    ax1.get_xaxis().set_tick_params(which='both', direction='out')
    ax1.get_yaxis().set_tick_params(which='both', direction='out')
    
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#006ddb' )
    
    for flier in bp['fliers']:
        flier.set(marker='o', color='#ff8000', alpha=0.5)
        
    for median in bp['medians']:
        median.set(color='#ff8000', linewidth=1)
    
    # color = '#006ddb'
    
    plt.xlabel('hamming distance')
    plt.ylabel('relative score')
    
    #ff8000
    #ax2 = plt.subplot(2,1,1)
    ax2 = plt.subplot(gs[0])
    plt.hist(distances, bins = np.arange(0, max(distances)+1.01,1), color = '#006ddb')
    ax2.get_xaxis().set_visible(False)
    
    ax2.get_xaxis().set_tick_params(bottom='off',labelbottom=False)
    ax2.get_yaxis().set_tick_params(which='both', direction='out')
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.locator_params(axis='y', nbins=4)
    plt.ylabel("frequency")
    plt.xlim([1, max(distances)+1])
    plt.tight_layout()
    
    plt.savefig(args.out_file + '.svg')
    plt.close()
    
    print('# Output file: %s'% args.out_file + '.svg')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Boxplot')
    parser.add_argument("-o", "--out_file", type=str, default= datetime.datetime.now().isoformat(), help='Name graphic')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.csv format. Separate multiple filenames with ";". If not set reads all *.csv in current directory')
    parser.add_argument("-d", "--directory", type = str, default=None, help='Set directory')
    args = parser.parse_args()

    box_plot(args)



