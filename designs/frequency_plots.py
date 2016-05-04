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
                                                                
    fig = plt.figure(figsize=(10,3))
            
    #ff8000
    ax2 = plt.subplot(1,1,1)
    plt.hist(distances, bins = np.arange(-0.5, max(distances)+1.5,1), color = '#006ddb')
    ax2.get_xaxis().set_tick_params(top='off', direction='out')
    ax2.get_yaxis().set_tick_params(which='both', direction='out')
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylabel("frequency")
    #plt.xlabel("hamming distance")
    plt.xlim([-0.5, max(distances)+0.5])
    plt.locator_params(axis='y', nbins=4)
    plt.locator_params(axis='x', nbins=32)
    
    plt.savefig(args.out_file + '.svg')
    plt.close()
    
    total = float((len(distances)))

    occurrences = dict((x, distances.count(x)) for x in set(distances))
    print("total " + str(total) )
    print("hamming distance" + "\t" +"count"+ "\t" + "percentage" )
    for entry in occurrences:
        count = float(occurrences[entry])
        percentage = count/total * 100
        print(str(entry) + "\t" + str(count) + "\t" + str(percentage))
      
    print('# Output file: %s'% args.out_file + '.svg')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Boxplot')
    parser.add_argument("-o", "--out_file", type=str, default= datetime.datetime.now().isoformat(), help='Name graphic')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.csv format. Separate multiple filenames with ";". If not set reads all *.csv in current directory')
    parser.add_argument("-d", "--directory", type = str, default=None, help='Set directory')
    args = parser.parse_args()

    box_plot(args)



