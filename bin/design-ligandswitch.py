#!/usr/bin/env python
# # -*- coding: utf-8 -*-

"""design-ligandswitch.py: Design a ligand triggered device."""

from __future__ import print_function

__author__ = "Stefan Hammer, Sven Findeiss"
__copyright__ = "Copyright 2018, Ribonets Project"
__credits__ = ["Stefan Hammer", "Sven Findeiss", "Ivo L. Hofacker",
                    "Christoph Flamm"]
__license__ = "GPLv3"
__version__ = "1.0"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hamer@univie.ac.at"
__status__ = "Production"

import RNAblueprint as rbp
import argparse
import sys
import time

try:
    from RNAsketch import *
except ImportError as e:
    print( "Error: %s" % e , file=sys.stderr)
    exit(1)


def main():
    parser = argparse.ArgumentParser(description='Design a ligand triggered device.')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-r", "--ratio", type=str, default='90:10', help='Ratio of the alternative to binding competent state in percent (default: 90:10)')
    parser.add_argument("-l", "--ligand", type=str, default="GAUACCAG&CCCUUGGCAGC;(...((((&)...)))...);-9.22", help='Binding motif and energy of the ligand (default: "GAUACCAG&CCCUUGGCAGC;(...((((&)...)))...);-9.22")')
    parser.add_argument("-T", "--temperature", type=float, default=37.0, help='Temperature of the energy calculations.')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    parser.add_argument("-s", "--stop", type=int, default=500, help='Stop optimization run if no better solution is aquired after (stop) trials.')
    parser.add_argument("-m", "--mode", type=str, default='random', help='Mode for getting a new sequence: sample, sample_plocal, sample_clocal, random')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    print("# Options: number={0:d}, stop={1:d}, mode={2:}, temperature={3:}, ratio={4:}, ligand={5:}".format(args.number, args.stop, args.mode, args.temperature, args.ratio, args.ligand))
    rbp.initialize_library(args.debug, args.kill)
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
        structures = ['(((((...((((((((.....)))))...)))...)))))........................',
            '.........................(((((((((((......)))))))))))...........']
        constraint = 'AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCANNNNNNNNNNNNNNNNNNNNNN'
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
    print("# " + "\n# ".join(structures) + "\n# " + constraint + "\n# " + start_sequence)

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
            print('# [' + str(i) + ']' + str(dg.component_vertices(i)) + ' -> ' + str(dg.number_of_sequences(i)) + ' solutions')

        # remember general DG values
        graph_properties = get_graph_properties(dg)
        # create a initial design object
        design = get_design(structures, start_sequence, constraint, args)

        # print header for csv file
        if (args.csv):
            print(";".join(["stop",
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
            design = get_design(structures, start_sequence, constraint, args)

            start = time.clock()
            try:
                (score, number_of_mutations) = adaptive_walk_optimization(dg, design, objective_function=ligand_objective, stop=args.stop, mode=args.mode, progress=args.progress)
            except ValueError as e:
                print (e.value)
                exit(1)

            ligand_objective(design, printDetails=True)
            # stop time counter
            sample_time = time.clock() - start

            if (args.csv):
                print(args.stop,
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

def get_design(structures, start_sequence, constraint, args):
    # split ligand motif
    ligandmotif = []
    try:
        ligandmotif = args.ligand.split(";")
        ligandmotif[-1] = float(ligandmotif[-1])
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)
        exit(1)

    design = vrnaDesign([], start_sequence)

    design.newState("pf", '.'*len(structures[0]))

    design.newState("bc", structures[0])
    design.state["bc"].constraint = structures[0]

    design.newState("ac", structures[1])
    design.state["ac"].constraint = structures[1]

    design.newState("pfl", '.'*len(structures[0]))
    design.state["pfl"].ligand=ligandmotif

    design.newState("bcl", structures[0])
    design.state["bcl"].constraint = structures[0]
    design.state["bcl"].ligand=ligandmotif

    # Set the given temperature for all states
    for state in design.state.values():
        state.temperature = args.temperature

    #add variables
    # split ratio
    design.ratio = []
    try:
        design.ratio = map(lambda x: (int(x) / float(100)), args.ratio.split(":")) #convert ratio to array of integers
        if sum(design.ratio) > 100:
            raise(ValueError("You specified ratios that in total become larger than 100 percent"))
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)
        exit(1)

    return design

def ligand_objective(design, printDetails=False):
    '''Calculates the objective function given a design object containing
    the partition function energies of the following four states.
    param design: RNAsektch design object with an additional design.ratio array
    returns: float score
    '''
    prob_bc_ligand = Z_from_G(design.state["bcl"].pf_energy - design.state["pfl"].pf_energy, design.state["pfl"].temperature)

    prob_ac = Z_from_G(design.state["ac"].pf_energy - design.state["pf"].pf_energy, design.state["pf"].temperature)

    prob_bc = Z_from_G(design.state["bc"].pf_energy - design.state["pf"].pf_energy, design.state["pf"].temperature)

    score = 1 - (prob_bc_ligand * (1 - abs(design.ratio[0] - prob_ac)) * (1 - abs(design.ratio[1] - prob_bc)))

    if printDetails:
        print("prob_bc_ligand: {0:}\nprob_ac: {1:}\nprob_bc: {2:}\nscore: {3:}".format(prob_bc_ligand, prob_ac, prob_bc, score))
    
    # we need to maximize the score
    return 1-score

def G_from_Z(Z, temperature):
    return - ((temperature + 273.15)*1.98717)/1000.0 * math.log(Z)

def Z_from_G(G, temperature):
    return math.exp(- (G/ (((temperature + 273.15)*1.98717)/1000.0)))


if __name__ == "__main__":
    main()
