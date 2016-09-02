from __future__ import print_function

try:
    from PyDesign import *
except ImportError, e:
    print(e.message)
    exit(1)

import RNAblueprint as rbp
import argparse
import sys
import time
import re

def main():
    parser = argparse.ArgumentParser(description='Design a multi-stable riboswitch similar to Hoehner 2013 paper.')
    parser.add_argument("-q", "--nupack", default=False, action='store_true', help='Use Nupack instead of the ViennaRNA package (for pseudoknots)')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    parser.add_argument("-e", "--exit", type=int, default=500, help='Exit optimization run if no better solution is aquired after (exit) trials.')
    parser.add_argument("-m", "--mode", type=str, default='random', help='Mode for getting a new sequence: sample, sample_local, sample_global, random')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    print("# Options: number={0:d}, exit={1:d}, mode={2:}, nupack={3:}".format(args.number, args.exit, args.mode, str(args.nupack)))
    rbp.initialize_library(args.debug, args.kill)
    # define structures
    structures = []
    temperatures = []
    constraint = ''
    start_sequence = ''
    
    data = ''
    for line in sys.stdin:
        data = data + '\n' + line
    
    if data:
        structures, constraint, start_sequence, temperatures = read_input_additions(data)
        temperatures = [float(t) for t in temperatures]
    else:
        structures = ['((((....))))....((((....))))........',
            '........((((....((((....))))....))))',
            '((((((((....))))((((....))))....))))']
        constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
        temperatures = [37.0, 37.0, 24.0]
    # try to construct dependency graph, catch errors and timeouts
    dg = None
    construction_time = 0.0
    sample_time = 0.0
        
    # construct dependency graph with these structures
    try:
        start = time.clock()
        dg = rbp.DependencyGraphMT(structures, constraint)
        construction_time = time.clock() - start
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)
    
    # general DG values
    print("# " + "\n# ".join(structures) + "\n# " + constraint)
    print("# Temperatures: ", temperatures)

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
            print(";".join(["exit",
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
            
            for i, t in enumerate(temperatures): 
                design.state[str(i)].temperature = t
                
                for j, s in enumerate(design.structures):
                    if j != i:
                        design.newState(str(j) + ':' + str(t), s, temperature=t)
            
            start = time.clock()
            
            # now do the optimization based on the chose mode for args.exit iterations
            try:
                (score, number_of_mutations) = classic_optimization(dg, design, objective_function=temp_objective, exit=args.exit, mode=args.mode, progress=args.progress)
            except ValueError as e:
                print (e.value)
                exit(1)
            # stop time counter
            sample_time = time.clock() - start
            
            if (args.csv):
                print(args.exit,
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

def temp_objective(design, weight=0.5):
    return calculate_objective_1(design) + weight * temp_objective_2(design)

def temp_objective_2(design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    (eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2) * 2 / (number_of_structures * (number_of_structures-1))
    
    :param design: Design object containing the sequence and structures
    :return: score calculated by the objective function
    '''
    objective_difference_part = 0
    
    for k in design.state.keys():
        # first iterate over all desired structures with their target temperature
        if re.match(re.compile("^[^\:]$"), k, flags=0):
            for kk in design.state.keys():
                if re.match(re.compile("^[^" + k + "]\:" + str(design.state[k].temperature)), kk, flags=0):
                    objective_difference_part += (design.state[k].eos - design.state[kk].eos)
                    print (" + ( " + k + " - " + kk + ")")
    print ("\n\n")
    if design.number_of_structures == 1:
        return objective_difference_part
    else:
        return objective_difference_part * 2 / (design.number_of_structures * (design.number_of_structures-1))

if __name__ == "__main__":
    main()


