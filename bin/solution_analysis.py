from __future__ import print_function

from collections import Counter
import argparse
import sys
import time
import re
import os
import gzip

import numpy as np
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description='Benchmark to check if the sampling is fair.')
    parser.add_argument("-f", "--file", type = str, default="", help='Read sequences from gzip file')
    parser.add_argument("-u", "--unfairfile", type = str, default="", help='Read unfair sequences from gzip file')
    parser.add_argument("-m", "--mode", type = str, default="fair", help='Set the mode to fair or unfair or both')
    args = parser.parse_args()

    print(";".join([
            "mode",
            "seq_length",
            "min_solution_count",
            "max_solution_count",
            "mean_solution_count",
            "number_found_solutions"]))
    
    unfair_values = []
    values = []
    
    with gzip.open(args.file, 'rb') as f:
        solutions = []
        for line in f:
            if line.startswith("#"):
                continue;
            else:
                solutions.append(line)
        
        solution_dict = Counter(solutions)
        values = solution_dict.values();
        mode = args.mode
        if args.mode == 'both':
            mode = 'fair'
        
        print(mode,
                    len(solutions[0]),
                    min(values),
                    max(values),
                    sum(values)/ float(len(values)), 
                    len(solution_dict.keys()), sep=";")
    
    # if it is both, make a plot with both curves
    if args.mode == 'both':
        with gzip.open(args.unfairfile, 'rb') as f:
            unfair_solutions = []
            for line in f:
                if line.startswith("#"):
                    continue;
                else:
                    unfair_solutions.append(line)
            
            unfair_solution_dict = Counter(unfair_solutions)
            unfair_values = unfair_solution_dict.values();
            
            print('unfair',
                        len(solutions[0]),
                        min(unfair_values),
                        max(unfair_values),
                        sum(unfair_values)/ float(len(unfair_values)), 
                        len(unfair_solution_dict.keys()), sep=";")
    
    
    plot_data(values, unfair_values, args);

def plot_data(data, unfair_data, args):
    fig = plt.figure()
    fig.set_size_inches(14, 7)
    plt.rc('font', family='Liberation Sans', size='18')
    if len(unfair_data) != 0:
        plt.hist(unfair_data, bins=max(unfair_data)-min(unfair_data), histtype='stepfilled', color='#FF8000', alpha=0.6, label='Unfair Sampling')
    plt.hist(data, bins=max(data)-min(data), histtype='stepfilled', color='#006DDB', alpha=0.7, label='Fair Sampling')
    plt.xlabel('Frequency of the solution found')
    plt.ylabel('Frequency')
    plt.yscale('symlog')
    plt.xscale('symlog')
    #plt.set_title('Histogram of Solution Counts')
    plt.grid(True)
    plt.legend(fontsize=18)
    plt.savefig(os.path.basename(args.file) + '.hist.svg')
    plt.clf()
    plt.close()
if __name__ == "__main__":
    main()

