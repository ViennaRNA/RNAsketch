from __future__ import print_function

from collections import Counter
import argparse
import sys
import time
import re
import os

def main():
    parser = argparse.ArgumentParser(description='Benchmark to check if the sampling is fair.')
    parser.add_argument("-f", "--file", type = str, default="", help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-u", "--unfair", default=False, action='store_true', help='Choose if you want to sample with the unfair libaray')
    parser.add_argument("-n", "--number", type=int, default=1000000, help='Number of times to sample a new solution')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    print("# Options: number={0:d}, kill={1:d}".format(args.number, args.kill))
    
    if (args.unfair):
        import RNAunfairdesign as rbp
    else:
        import RNAblueprint as rbp
    
    # initialize library
    rbp.initialize_library(args.debug, args.kill)
    
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
    
    # try to construct dependency graph, catch errors and timeouts
    dg = None
    
    try:
        dg = rbp.DependencyGraphMT(structures, constraint)
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)
    
    # general DG values
    if (dg is not None):
        print("# " + "\n# ".join(structures) + "\n# " + constraint)
        # print the amount of solutions
        print('# Maximal number of solutions: ' + str(dg.number_of_sequences()))
        # print the amount of connected components
        number_of_components = dg.number_of_connected_components()
        print('# Number of Connected Components: ' + str(number_of_components))
        for i in range(0, number_of_components):
            print('# [' + str(i) + ']' + str(dg.component_vertices(i)))
        
        # if requested write out a graphml file
        if args.graphml is not None:
            with open(args.graphml, 'w') as f:
                f.write(dg.get_graphml() + "\n")

        # main purpose, count the sequences
        for n in range(0, args.number):
            dg.sample()
            print (dg.get_sequence())
            if (args.progress):
                sys.stderr.write("\r{0:5.2f}%".format(float(n)/float(args.number)*100))
                sys.stderr.flush()
        if (args.progress):
            sys.stdout.write("\r" + " " * 60 + "\r")
            sys.stdout.flush()

if __name__ == "__main__":
    main()

