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
        for i in range(0, number_of_components):
            print('# [' + str(i) + ']' + str(dg.component_vertices(i)))

        # remember general DG values
        graph_properties = get_graph_properties(dg)
        # create a initial design object

        if (args.nupack):
            design = nupackDesign(structures, start_sequence)
        else:
            design = vrnaDesign(structures, start_sequence)

        if not design.sequence:
            dg.sample()
            design.sequence = dg.get_sequence()
        else:
            dg.set_sequence(design.sequence)


        # calculations for plot

        # calculate objectives for initial sequence
        x_center = calculate_objective_1(design)
        y_center = calculate_objective_2(design)

        x_origin = 0 # = x_center-x_center
        y_origin = 0

        plt.xlabel('objective 2')
        plt.ylabel('objective 1')
        plt.plot([x_origin], [y_origin], 'ro')

        ax_lim = 0

        for i in range(0, args.number):
            PyDesign._sample_sequence(dg, design, args.mode, args.sample_steps)
            x_new = calculate_objective_2(design)
            y_new = calculate_objective_1(design)

            # calculate relative positions to initial sequence
            x = x_new - x_center
            y = y_new - y_center

            # find largest x/y and use as limit for both axes
            if abs(x) > ax_lim:
                ax_lim = abs(x)
            if abs(y) > ax_lim:
                ax_lim = abs(y)

            plt.plot([x], [y], 'bo')

        plt.axhline(color='k')
        plt.axvline(color='k')

        ax_lim += 1 # +1 to avoid points on axes
        plt.xlim([-ax_lim, ax_lim])
        plt.ylim([-ax_lim, ax_lim])

        plt.savefig(args.out_file)
        print('# Output file: %s'% args.out_file)
        plt.close()

    else:
        print('# Construction time out reached!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot...')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-m", "--mode", type=str, default='sample_global', help='Mode for getting a new sequence: sample, sample_local, sample_global, sample_strelem')
    parser.add_argument("-s", "--sample_steps", type=int, default=1, help='Count how many times to do the sample operation')
    parser.add_argument("-n", "--number", type=int, default=100, help='Count how many times to sample sequence')
    parser.add_argument("-q", "--nupack", default=False, action='store_true', help='Use Nupack instead of the ViennaRNA package (for pseudoknots)')
    parser.add_argument("-o", "--out_file", type=str, default= datetime.datetime.now().isoformat()+".png", help='Name output file')
    args = parser.parse_args()

    plot_sequence_objective(args)



