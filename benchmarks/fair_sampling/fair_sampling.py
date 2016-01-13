from __future__ import print_function

from collections import Counter
import argparse
import sys
import time
import re
import os

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def main():
    parser = argparse.ArgumentParser(description='Benchmark to check if the sampling is fair.')
    parser.add_argument("-f", "--file", type = str, default="", help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-n", "--number", type=int, default=1000000, help='Number of times to sample a new solution')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    print("# Options: number={0:d}, kill={1:d}".format(args.number, args.kill))
    
    # define structures
    structures = []
    constraint = ""
    if (args.input):
        for line in sys.stdin:
            if re.match(re.compile("[\(\)\.]"), line, flags=0):
                structures.append(line.rstrip('\n'))
            elif re.match(re.compile("[ACGTUWSMKRYBDHVN]"), line, flags=0):
                constraint = line.rstrip('\n')
            elif re.search(re.compile("@"), line, flags=0):
                break;
    elif (args.file is not ""):
        print("# Input File: {0:}".format(os.path.basename(args.file)))
        with open(args.file) as f:
            data = f.read()
            lines = data.split("\n")
            for line in lines:
                if re.match(re.compile("[\(\)\.]"), line):
                    structures.append(line)
                if re.match(re.compile("[\ AUGC]"), line):
                    elements = str(line)
                    constraint = elements.replace(" ", "N")
                if line.startswith(";"):
                    break
    else:
        structures = ['((((....))))....((((....))))........',
            '........((((....((((....))))....))))',
            '((((((((....))))((((....))))....))))']
        constraint = 'AGNYNCNNNGNWNANGCNCNNCAUNNYNRNANAGNN'

    print(";".join([
            "mode",
            "seq_length",
            "graph_construction",
            "num_cc",
            "max_specials",
            "max_component_vertices",
            "max_special_ratio",
            "mean_special_ratio",
            "nos",
            "construction_time",
            "sample_time",
            "min_solution_count",
            "max_solution_count",
            "mean_solution_count",
            "number_found_solutions"]))
    # store the values
    unfair_values = []
    values = []
       
    # do everything twice with fair and unfair sampling
    import RNAdesign as rd
    for mode in ['fair', 'unfair']:
        rd.initialize_library(args.debug, args.kill)
        
        # try to construct dependency graph, catch errors and timeouts
        dg = None
        graph_construction = 0
        construction_time = 0.0
        sample_time = 0.0
        max_specials = 0
        max_component_vertices = 0
        max_special_ratio = 0
        mean_special_ratio = 0
        num_cc = 0
        nos = 0
        
        try:
            start = time.clock()
            dg = rd.DependencyGraphMT(structures, constraint)
            construction_time = time.clock() - start
        except Exception as e:
            print( "Error: %s" % e , file=sys.stderr)
        
        # general DG values
        if (dg is not None):
            if mode is 'fair':
                print("# " + "\n# ".join(structures) + "\n# " + constraint)
                # print the amount of solutions
                print('# Maximal number of solutions: ' + str(dg.number_of_sequences()))
                # print the amount of connected components
                number_of_components = dg.number_of_connected_components()
                print('# Number of Connected Components: ' + str(number_of_components))
                for i in range(0, number_of_components):
                    print('# [' + str(i) + ']' + str(dg.component_vertices(i)))
            
            # remember general DG values
            graph_construction = 1
            num_cc = dg.number_of_connected_components()
            nos = dg.number_of_sequences()
            special_ratios = []
            for cc in range(0, num_cc):
                cv = len(dg.component_vertices(cc))
                sv = len(dg.special_vertices(cc))
                special_ratios.append(float(sv)/float(cv))
                if (max_specials < sv):
                    max_specials = sv
                if (max_component_vertices < cv):
                    max_component_vertices = cv
            max_special_ratio = max(special_ratios)
            mean_special_ratio = sum(special_ratios)/len(special_ratios)

            # if requested write out a graphml file
            if args.graphml is not None:
                with open(args.graphml, 'w') as f:
                    f.write(dg.get_graphml() + "\n")

            # main purpose, count the sequences
            start = time.clock()
            result = count_sequences(dg, mode, args)
            sample_time = time.clock() - start
            
            if mode == 'fair':
                values = result[4]
            else:
                unfair_values = result[4]
            
            print(mode,
                        len(structures[0]),
                        graph_construction,
                        num_cc,
                        max_specials,
                        max_component_vertices,
                        max_special_ratio,
                        mean_special_ratio,
                        nos,
                        construction_time,
                        sample_time, 
                        result[0],
                        result[1],
                        result[2], 
                        result[3], sep=";")
            # now once more with unfair sampling
            import RNAunfairdesign as rd
    
    plot_data(values, unfair_values, args);

def count_sequences(dg, mode, args):
    solutions = []
    for n in range(0, args.number):
        dg.sample()
        solutions.append(dg.get_sequence())
        if (args.progress):
            sys.stdout.write("\r{0:}: {1:5.2f}%".format(mode, float(n)/float(args.number)*100))
            sys.stdout.flush()
    if (args.progress):
        sys.stdout.write("\r" + " " * 60 + "\r")
        sys.stdout.flush()
    solution_dict = Counter(solutions)
    if args.debug:
        for k,v in solution_dict.items():
            print('#', k, ": ", v)
    values = solution_dict.values();
    return [ min(values), max(values), sum(values)/ float(len(values)), len(solution_dict.keys()), values ]

def plot_data(data, unfair_data, args):
    #plt.hist(data, normed=True, facecolor='lightblue', bins=max(data)-min(data), histtype='stepfilled', label='fair')
    if len(unfair_data) != 0:
        plt.hist(unfair_data, bins=max(unfair_data)-min(unfair_data), histtype='stepfilled', color='red', alpha=0.6, label='unfair')
    plt.hist(data, bins=max(data)-min(data), histtype='stepfilled', color='lightgreen', alpha=0.7, label='fair')
    plt.xlabel('Solution Count')
    plt.ylabel('Frequency')
    plt.yscale('symlog')
    plt.xscale('symlog')
    plt.title('Sampling: ' + os.path.basename(args.file))
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.basename(args.file) + '.hist.svg')
    plt.clf()
    plt.close()

if __name__ == "__main__":
    main()

