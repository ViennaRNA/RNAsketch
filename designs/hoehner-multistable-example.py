from __future__ import print_function

from PyDesign import *
import RNAdesign as rd
import argparse
import sys
import time

def main():
    parser = argparse.ArgumentParser(description='Design a tri-stable example same to Hoehner 2013 paper.')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    parser.add_argument("-j", "--jump", type=int, default=1000, help='Do random jumps in the solution space for the first (jump) trials.')
    parser.add_argument("-e", "--exit", type=int, default=1000, help='Exit optimization run if no better solution is aquired after (exit) trials.')
    parser.add_argument("-m", "--mode", type=str, default='sample_global', help='Mode for getting a new sequence: sample, sample_local, sample_global')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    print("# Options: number={0:d}, jump={1:d}, exit={2:d}, mode={3:}".format(args.number, args.jump, args.exit, args.mode))
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
    
    # create the design object
    design = Design(structures, start_sequence)
    
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
    
    # print header for csv file
    if (args.csv):
        print(";".join(["jump",
                    "exit",
                    "mode",
                    "score",
                    "num_mutations", 
                    "graph_construction",
                    "num_cc",
                    "max_specials",
                    "max_component_vertices",
                    "max_special_ratio",
                    "mean_special_ratio",
                    "nos",
                    "construction_time",
                    "sample_time",
                    design.write_csv_header()]))
        
    # construct dependency graph with these structures
    try:
        start = time.clock()
        dg = rd.DependencyGraphMT(structures, constraint)
        construction_time = time.clock() - start
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)

    if (dg is not None):
        # general DG values
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

        # main loop from zero to number of solutions
        for n in range(0, args.number):
            start = time.clock()
            (score, number_of_mutations) = optimization_run(dg, design, args)
            sample_time = time.clock() - start
            
            if (args.progress):
                sys.stdout.write("\r" + " " * 60 + "\r")
                sys.stdout.flush()
            
            if (args.csv):
                print(args.jump,
                        args.exit,
                        "\"" + args.mode + "\"",
                        score,
                        number_of_mutations,
                        graph_construction,
                        num_cc,
                        max_specials,
                        max_component_vertices,
                        max_special_ratio,
                        mean_special_ratio,
                        nos,
                        construction_time,
                        sample_time,
                        design.write_csv(), sep=";")
            else:
                print(design.write_out(score))
    else:
        print(args.jump,
                args.exit,
                "\"" + args.mode + "\"",
                0,
                0,
                graph_construction,
                num_cc,
                max_specials,
                max_component_vertices,
                max_special_ratio,
                mean_special_ratio,
                nos,
                construction_time,
                sample_time, sep=";")

# main optimization
def optimization_run(dg, design, args):
    score = 0
    count = 0
    jumps = args.jump
    # sample scratch sequence
    dg.sample()
    # print this sequence with score
    design.sequence = dg.get_sequence()
    score = calculate_objective(design);
    #print dg.get_sequence() + '\t' + str(score)
    
    # sample globally for num_opt times and print
    number_of_mutations = 0
    while 1:
        # sample sequence
        if jumps:
            mut_nos = dg.sample()
            jumps -= 1
            count = 0
        else:
            if args.mode == 'sample':
                mut_nos = dg.sample()
            elif args.mode == 'sample_global':
                mut_nos = dg.sample_global()
            elif args.mode == 'sample_local':
                mut_nos = dg.sample_local()
            else:
                sys.stdout.write("Wrong sample argument: " + args.mode + "\n")
                sys.exit(1)
        # write progress
        if (args.progress):
            sys.stdout.write("\rMutate: {0:7.0f}/{1:5.0f} from NOS: {2:7.0f}".format(number_of_mutations, count, mut_nos) + " " * 20)
            sys.stdout.flush()
        
        design.sequence = dg.get_sequence()
        this_score = calculate_objective(design);
        
        if (this_score < score):
            score = this_score
            count = 0
        else:
            dg.revert_sequence();
            design.sequence = dg.get_sequence()
            count += 1
            if count > args.exit:
                break
        number_of_mutations += 1
    
    # finally return the result
    return score, number_of_mutations

if __name__ == "__main__":
    main()


