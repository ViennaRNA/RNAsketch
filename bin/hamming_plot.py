from __future__ import print_function

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
import math

def plot_hamming(args):
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
    
    distances = []
    for filename in filenames:
        with open(filename) as sampled_seq:
            reader = csv.reader(sampled_seq, delimiter=";")
            header = reader.next()

            # get information of intial sequence
            initial_sequence = reader.next()
            for line in reader:
            
                distances.append(int(line[3]))
                
        print(distances)
        plt.hist(distances, bins = np.arange(-0.5, max(distances)+1.5,1), color = '#006ddb')
        plt.title("Hamming Distance to Neighbors")
        plt.xlabel("hamming distance")
        plt.ylabel("frequency")

        fig = plt.gcf()

        plt.savefig(filename + '.svg')
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot...')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.csv format. Separate multiple filenames with ";". If not set reads all *.csv in current directory')
    parser.add_argument("-d", "--directory", type = str, default=None, help='Set directory')
    args = parser.parse_args()

    plot_hamming(args)
