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
import math


def plot_sequence_objective(args):

    # define structures
    structures = []
    constraint = ''
    start_sequence = ''

    if (args.input):
        data = ''
        for line in sys.stdin:
            data = data + '\n' + line
        (structures, constraint, start_sequence) = read_input(data)
    elif (args.file is not None):
        print("# Input File: {0:}".format(args.file))
        (structures, constraint, start_sequence) = read_inp_file(args.file)
    else:
        structures = ['((((....))))....((((....))))........',
            '........((((....((((....))))....))))',
            '((((((((....))))((((....))))....))))']
        constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    # try to construct dependency graph, catch errors and timeouts
    dg = None
    construction_time = 0.0
    sample_time = 0.0


    # construct dependency graph with these structures
    try:
        start = time.clock()
        dg = rd.DependencyGraphMT(structures, constraint)
        construction_time = time.clock() - start
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)

    # general DG values
    print("# " + "\n# ".join(structures) + "\n# " + constraint)

    if (dg is not None):

        # print the amount of solutions
        print('# Maximal number of solutions: ' + str(dg.number_of_sequences()))
        # print the amount of connected components
        number_of_components = dg.number_of_connected_components()
        print('# Number of Connected Components: ' + str(number_of_components))
        lnos_sum = 0
        for i in range(0, number_of_components):
            lnos = dg.sample_global(i)
            lnos_sum += lnos
            dg.revert_sequence()
            print('# [' + str(i) + ']' + str(dg.component_vertices(i)) + '\t' + str(lnos))
           
        print('LNOS: ' + str(lnos_sum))
        
        number = args.number
        perc_LNOS = args.percentage_of_LNOS
        
        print('LNOS * percentage: ' + str(int(math.ceil(lnos_sum * perc_LNOS))))

        print('number: ' + str(number))
        if number > lnos_sum * perc_LNOS:
            print('Number > 0.85 * LNOS!')
            #number = int(math.ceil(lnos_sum * perc_LNOS))
            
        if number == -1:
            number = int(math.ceil(lnos_sum * perc_LNOS))
                
        # remember general DG values
        graph_properties = get_graph_properties(dg)

        start_sequences = [] 
        
        if args.start_seq is not None:
            
            with open(args.start_seq) as f:
                data = f.read()
            start_seq = data.split("\n")  
            for seq in start_seq:
                if seq != '':
                    start_sequences.append(seq.rstrip('\n'))
       
        if args.start_seq is None:
            start_sequences.append(start_sequence)
        
        for i in range(0,len(start_sequences)):  
        
            if args.start_seq is not None:   
                csv_file_name = args.out_file + "_" + args.mode + "_start_seq_" + str(i + 1)   
            else:
                csv_file_name = args.out_file + "_" + args.mode  
                
            csv_file = open(csv_file_name + ".csv", 'w') 
            start_sequence = start_sequences[i]

            # create an initial design object
            if (args.nupack):
                design = nupackDesign(structures, start_sequence)
            else:
                design = vrnaDesign(structures, start_sequence)

            if not design.sequence:
                dg.sample()
                design.sequence = dg.get_sequence()
            else:
                dg.set_sequence(design.sequence)

            # print header for csv file                           
            csv_file.write(";".join(["objective 1",
                            "objective 2",
                            "score",
                            "hamming distance",
                            "sequence \n"]))

            # calculate objectives for initial sequence
            x_center = calculate_objective_1(design)
            y_center = calculate_objective_2(design)

            # print csv input for initial sequence
            score = x_center + args.weight * y_center

            # store entries -> just unique entries
            samples = set()
            samples.add(design.sequence)
            
            initial_seq = design.sequence
            csv_file.write(";".join([str(x_center),
                    str(y_center),
                    str(score),
                    str(RNA.hamming_distance(initial_seq, initial_seq)),
                    design.sequence]) + "\n")
                    
            for i in range(0, number): 

                if args.exit is not None:
                    no_tries = 0
                    
                while design.sequence in samples:
                    (mut_nos, sample_count) = PyDesign.sample_sequence(dg, design, args.mode, args.sample_steps)
                    dg.revert_sequence(sample_count)
                    
                    if args.exit is not None:
                        no_tries += 1

                        if no_tries == args.exit:
                            break
                
                if args.exit is not None and no_tries == args.exit:
                    break
                    
                x_new = calculate_objective_1(design)
                y_new = calculate_objective_2(design)


                if args.progress:
                    sys.stdout.write("\r# Sampling: {0:7.0f}/{1:5.0f}".format(i + 1, number) + " " * 20)
                    sys.stdout.flush()

                score = x_new + args.weight * y_new


                csv_file.write(";".join([str(x_new),
                    str(y_new),
                    str(score),
                    str(RNA.hamming_distance(initial_seq, design.sequence)),
                    design.sequence]) + "\n")

                samples.add(design.sequence)
                
            print('\n # Output file: %s'% csv_file_name + '.csv')

    else:
            print('# Construction time out reached!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Uses PyDesign.sample_sequence to explore sequence space and writes obj1,obj2,score and hamming distance (relative to first sequence) in *.csv')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-y", "--start_seq", type = str, default=None, help='Read start_seq')
    parser.add_argument("-m", "--mode", type=str, default='sample_global', help='Mode for getting a new sequence: sample, sample_local, sample_global, sample_strelem')
    parser.add_argument("-s", "--sample_steps", type=int, default=1, help='Count how many times to do the sample operation')
    parser.add_argument("-n", "--number", type=int, default=1000, help='Define number of neighbors')
    parser.add_argument("-e", "--exit", type=int, default=None, help='Exit value')
    parser.add_argument("-l", "--percentage_of_LNOS", type=float, default=0.85, help='Define percentage of LNOS')
    parser.add_argument("-w", "--weight", type=float, default=0.5, help='Define weighting-factor')
    parser.add_argument("-q", "--nupack", default=False, action='store_true', help='Use Nupack instead of the ViennaRNA package (for pseudoknots)')
    parser.add_argument("-o", "--out_file", type=str, default= datetime.datetime.now().isoformat(), help='Name file')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    args = parser.parse_args()

    plot_sequence_objective(args)



