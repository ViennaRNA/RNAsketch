from __future__ import print_function

from collections import Counter
import argparse
import sys
import time
import re
import os
import gzip

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def main():
    parser = argparse.ArgumentParser(description='Benchmark to check if the sampling is fair.')
    parser.add_argument("-f", "--file", type = str, default="", help='Read sequences from gzip file')
    parser.add_argument("-m", "--mode", type = str, default="fair", help='Set the mode to fair or unfair')
    args = parser.parse_args()

    print(";".join([
            "mode",
            "seq_length",
            "min_solution_count",
            "max_solution_count",
            "mean_solution_count",
            "number_found_solutions"]))
    
    with gzip.open(args.file, 'rb') as f:
        solutions = []
        for line in f:
            if line.startswith("#"):
                continue;
            else:
                solutions.append(line)
        
        solution_dict = Counter(solutions)
        values = solution_dict.values();
        plot_data(values, args);
        
        print(args.mode,
                    len(solutions[0]),
                    min(values),
                    max(values),
                    sum(values)/ float(len(values)), 
                    len(solution_dict.keys()), sep=";")

def plot_data(data, args):
    plt.hist(data, normed=True, facecolor='lightblue', bins=60, histtype='step')
    plt.xlabel('Solution Count')
    plt.ylabel('Frequency')
    plt.title(args.mode + ' sampling: ' + os.path.basename(args.file))
    plt.grid(True)
    plt.savefig(os.path.basename(args.file) + '.' + args.mode + '.svg')
    plt.clf()
    plt.close()

if __name__ == "__main__":
    main()

