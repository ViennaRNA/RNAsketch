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

def main():
    parser = argparse.ArgumentParser(description='Design a tri-stable example same to Hoehner 2013 paper.')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-q", "--nupack", default=False, action='store_true', help='Use Nupack instead of the ViennaRNA package (for pseudoknots)')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    parser.add_argument("-j", "--jump", type=int, default=300, help='Do random jumps in the solution space for the first (jump) trials.')
    parser.add_argument("-e", "--exit", type=int, default=500, help='Exit optimization run if no better solution is aquired after (exit) trials.')
    parser.add_argument("-s", "--strelem", type=int, default=1800, help='Optimize structural elements and exit after (strelem) unsucessful trials.')
    parser.add_argument("-m", "--mode", type=str, default='sample_global', help='Mode for getting a new sequence: sample, sample_local, sample_global, sample_strelem')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    print("# Options: number={0:d}, jump={1:d}, exit={2:d}, strelem={3:d}, mode={4:}, nupack={5:}".format(args.number, args.jump, args.exit, args.strelem, args.mode, str(args.nupack)))
    rd.initialize_library(args.debug, args.kill)
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
        
        # if requested write out a graphml file
        if args.graphml is not None:
            with open(args.graphml, 'w') as f:
                f.write(dg.get_graphml() + "\n")
        
        # print the amount of solutions
        print('# Maximal number of solutions: ' + str(dg.number_of_sequences()))
        # print the amount of connected components
        number_of_components = dg.number_of_connected_components()
        print('# Number of Connected Components: ' + str(number_of_components))
        for i in range(0, number_of_components):
            print('# [' + str(i) + ']' + str(dg.component_vertices(i)))
        
        # remember general DG values
        graph_properties = get_graph_properties(dg)
        # create a initial design object
        if (args.nupack):
            design = nupackDesign(structures, start_sequence)
        else:
            design = vrnaDesign(structures, start_sequence)
        
        # print header for csv file
        if (args.csv):
            print(";".join(["jump",
                        "exit",
                        "strelem",
                        "mode",
                        "score",
                        "num_mutations",
                        "construction_time",
                        "sample_time",
                        design.write_csv_header()] +
                        graph_properties.keys()))

        # main loop from zero to number of solutions
        for n in range(0, args.number):
            # reset the design object
            if (args.nupack):
                design = nupackDesign(structures, start_sequence)
            else:
                design = vrnaDesign(structures, start_sequence)
            
            start = time.clock()
            # do a complete sampling jump times
            #(score, number_of_jumps) = classic_optimization(dg, design, exit=args.jump, mode='sample', progress=args.progress)
            # now do the optimization based on the chose mode
            try:
                (score, number_of_mutations) = classic_optimization_const(dg, design, exit=args.exit, mode=args.mode, progress=args.progress)
            except ValueError as e:
                print (e.value)
                exit(1)
            # now do the optimization with mode strelem where we take structural elements and replace them a little
            number_of_strelem = 0
            #if forgi_available:
            #    (score, number_of_strelem) = classic_optimization(dg, design, exit=args.strelem, mode='sample_strelem', progress=args.progress)
            #else:
            #    sys.stderr.write("-" * 60 + "\nWARNING: Strelem sampling not available!!!\nPlease install forgi https://github.com/pkerpedjiev/forgi\n" + "-" * 60 + "\n")
            #    sys.stderr.flush() 
            # sum up for a complete number of mutations
            #number_of_mutations += number_of_jumps + number_of_strelem
            sample_time = time.clock() - start
            
            if (args.csv):
                print(args.jump,
                        args.exit,
                        args.strelem,
                        "\"" + args.mode + "\"",
                        score,
                        number_of_mutations,
                        construction_time,
                        sample_time,
                        design.write_csv(),
                        *graph_properties.values(), sep=";")
            else:
                print(design.write_out(score))
    else:
        print('# Construction time out reached!')


def classic_optimization_const(dg, design, objective_function=calculate_objective, exit=1000, mode='sample', progress=False):
    '''
    Takes a Design object and does a classic optimization of this sequence.
    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_functions: array of functions which takes a design object and returns a score for evaluation
    :param exit: Number of unsuccessful new sequences before exiting the optimization
    :param mode: String defining the sampling mode: sample, sample_global, sample_local
    :param progress: Whether or not to print the progress to the console
    :param return: Optimization score reached for the final sequence
    "param return: Number of samples neccessary to reach this result
    '''
    # if the design has no sequence yet, sample one from scratch
    if not design.sequence:
        dg.sample()
        design.sequence = dg.get_sequence()
    else:
        dg.set_sequence(design.sequence)
    
    score = objective_function(design)
    # count for exit condition
    count = 0
    # remember how may mutations were done
    number_of_samples = 0
    
    # main optimization loop 
    while True:
        # count up the mutations
        number_of_samples += 1
        # exit condition
        if number_of_samples > exit:
            break
        # sample a new sequence
        (mut_nos, sample_count) = PyDesign._sample_sequence(dg, design, mode)
        
        # write progress
        if progress:
            sys.stderr.write("\rMutate: {0:7.0f}/{1:5.0f} | Score: {2:5.2f} | NOS: {3:.5e}".format(number_of_samples, count, score, mut_nos) + " " * 20)
            sys.stderr.flush()
        
        this_score = objective_function(design)
        # evaluate
        if (this_score < score):
            score = this_score
            count = 0
        else:
            dg.revert_sequence(sample_count)
            design.sequence = dg.get_sequence()
            count += 1
    
    # clear the console
    if (progress):
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()
    
    # finally return the result
    return score, number_of_samples


if __name__ == "__main__":
    main()


