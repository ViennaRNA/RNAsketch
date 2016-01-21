from __future__ import print_function

import RNAdesign as rd
import RNA
import argparse
import sys
import collections
from collections import deque
import re
import math
import time

kT = ((37+273.15)*1.98717)/1000.0; # kT = (betaScale*((temperature+K0)*GASCONST))/1000.0; /* in Kcal */

class Result:
    def __init__(self, sequence, score, structures, number_of_mutations):
        self.sequence = sequence
        self.score = score
        self.structures = structures
        self.number_of_mutations = number_of_mutations
        self.eos = []
        self.probs = []
        
        (self.mfe_struct, self.mfe_energy) = RNA.fold(self.sequence)
        self.part_funct = RNA.pf_fold(self.sequence)[1]
        for struct in self.structures:
            this_eos = RNA.energy_of_struct(self.sequence, struct)
            self.eos.append(this_eos)
            self.probs.append( math.exp((self.part_funct-this_eos) / kT ) )
    def write_out(self):
        #first clean up last line
        sys.stdout.write("\r" + " " * 60 + "\r")
        sys.stdout.flush()
        print(self.sequence + '\t{0:9.4f}'.format(self.score))
        for i, struct in enumerate(self.structures):
            print(struct + '\t{0:9.4f}\t{1:+9.4f}\t{2:9.4f}'.format(self.eos[i], self.eos[i]-self.mfe_energy, self.probs[i]))
        print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))


def main():
    parser = argparse.ArgumentParser(description='constrained sequence generation')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    parser.add_argument("-s", "--size_constraint", type=int, default=100, help='Size of negative constraints container')
    parser.add_argument("-e", "--exit", type=int, default=1000, help='Exit optimization run if no better solution is aquired after exit trials.')
    parser.add_argument("-m", "--mode", type=str, default='sample_global', help='Mode for getting a new sequence: sample, sample_local, sample_global')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    parser.add_argument("-b", "--begin", type=str, default=None, help='set sequence to begin the optimization with')
    args = parser.parse_args()
    
    print("# Options: number={0:d}, size_constraint={1:d}, exit={2:d}, mode={3:}".format(args.number, args.size_constraint, args.exit, args.mode))
    rd.initialize_library(args.debug, args.kill)

    # define structural and sequence constraints
    pos_constraint = []
    seq_constraint = ""
    if (args.input):
        for line in sys.stdin:
            if re.match(re.compile("[\(\)\.]"), line, flags=0):
                pos_constraint.append(line.rstrip('\n'))
            elif re.match(re.compile("[ACGTUWSMKRYBDHVN]"), line, flags=0):
                seq_constraint = line.rstrip('\n')
            elif re.search(re.compile("@"), line, flags=0):
                break;
    elif (args.file is not None):
        print("# Input File: {0:}".format(args.file))
        with open(args.file) as f:
            data = f.read()
            lines = data.split("\n")
            for line in lines:
                if re.match(re.compile("[\(\)\.]"), line):
                    pos_constraint.append(line)
                if re.match(re.compile("[\ AUGC]"), line):
                    elements = str(line)
                    seq_constraint = elements.replace(" ", "N")
                if line.startswith(";"):
                    break
    else:
        pos_constraint = ['((((....))))....((((....))))........',
            '........((((....((((....))))....))))',
            '((((((((....))))((((....))))....))))']
        seq_constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    
    # generate negative constraints container
    neg_structures = collections.deque(maxlen=args.size_constraint)
    
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
        mfe_reached_str = ""
        diff_eos_mfe_str = ""
        prob_str = ""
        for s in range(0, len(pos_constraint)):
            mfe_reached_str = mfe_reached_str + "mfe_reached_" + str(s) +";"
            diff_eos_mfe_str = diff_eos_mfe_str + "diff_eos_mfe_" + str(s) + ";"
            prob_str = prob_str + "prob_" + str(s) + ";"
        print(";".join(["size_constraint",
                    "exit",
                    "mode",
                    "score",
                    "num_mutations", 
                    "seq_length",
                    "sequence",
                    "graph_construction",
                    "num_cc",
                    "max_specials",
                    "max_component_vertices",
                    "max_special_ratio",
                    "mean_special_ratio",
                    "nos",
                    "construction_time",
                    "sample_time"]) + ";" + 
                    mfe_reached_str + 
                    diff_eos_mfe_str +
                    prob_str)
    
    # construct dependency graph with the pos_constraint
    try:
        start = time.clock()
        dg = rd.DependencyGraphMT(pos_constraint, seq_constraint)
        construction_time = time.clock() - start
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)
    
    # main loop from zero to number of solutions
    if (dg is not None):
        # general DG values
        print("# " + "\n# ".join(pos_constraint) + "\n# " + seq_constraint)
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
        for i in range(0, args.number):
            start = time.clock()
            r = optimization(dg, pos_constraint, neg_structures, args)
            sample_time = time.clock() - start
            if (args.csv):
                # process result and write result of this optimization to stdout
                diff_eos_mfe = []
                mfe_reached = []
                for i in range(0, len(r.structures)):
                    mfe_reached.append(0)
                    eos_mfe = r.eos[i] - r.mfe_energy
                    diff_eos_mfe.append(eos_mfe)
                    if r.eos[i] == r.mfe_energy:
                        mfe_reached[i] = 1

                if (args.progress):
                    sys.stdout.write("\r" + " " * 60 + "\r")
                    sys.stdout.flush()

                print(args.size_constraint,
                        args.exit,
                        "\"" + args.mode + "\"",
                        r.score,
                        r.number_of_mutations, 
                        len(r.sequence),
                        "\"" + r.sequence + "\"",
                        graph_construction,
                        num_cc,
                        max_specials,
                        max_component_vertices,
                        max_special_ratio,
                        mean_special_ratio,
                        nos,
                        construction_time,
                        sample_time,
                        *(mfe_reached + diff_eos_mfe + r.probs), sep=";")
            else:
                r.write_out()
    else:
        dummylist = [0] * len(pos_constraint)
        print(args.size_constraint,
                args.exit,
                "\"" + args.mode + "\"",
                0, 
                0,
                len(pos_constraint[0]),
                "\"\"",
                graph_construction,
                num_cc,
                max_specials,
                max_component_vertices,
                max_special_ratio,
                mean_special_ratio,
                nos,
                construction_time,
                sample_time,
                *(dummylist + dummylist + dummylist), sep=";")


def optimization(dg, pos_constraint, neg_structures, args):
    score = float('Infinity')
    count = 0
    
    # get the initial sequence
    current_seq = dg.get_sequence()
    if(args.begin is not None)
        dg.set_sequence("")
    else
        dg.sample()        
    # number of mutations
    i = 0

    while True:
        (current_seq_struct, current_seq_mfe) = RNA.fold(current_seq)
        # if shapes option on -> calculate shape:
        # RNAshapes -D '(((((....)))))..((...((....))))....' -t 1
        # and then string-compare shapes 
        
        # if tree-distance optin is on calculate tree-edit distance between
        # current struct and all pos constraint. if dist < args.treexxx then ok 
        # int = RNA.tree_edit_distance(RNA.make_tree(struct1), RNA.make_tree(struct2))
        if current_seq_struct in pos_constraint:

            current_score = calculate_objective(current_seq, pos_constraint)

            if current_score < score:
                score = current_score
                result_seq = current_seq
                count = 0

            else:
                count += 1
                if count > args.exit:
                    break

        if current_seq_struct not in neg_structures:
            neg_structures.append(current_seq_struct)

        while True:
            if args.mode == 'sample':
                mut_nos = dg.sample()
            elif args.mode == 'sample_global':
                mut_nos = dg.sample_global()
            elif args.mode == 'sample_local':
                mut_nos = dg.sample_local()
            else:
                sys.stdout.write("Wrong sample argument: " + args.mode + "\n")
                sys.exit(1)
            current_seq = dg.get_sequence()
            perfect = True

            if args.progress:
                sys.stdout.write("\rMutate: {0:5.0f}/{1:5.0f} from NOS:{2:7.0f}".format(count, i, mut_nos))
                sys.stdout.flush()
            #calculate eos of neg and pos constraints before, because we need it again and again
            neg_energy_constraint = []
            for n in range(0, len(neg_structures)):
                neg_energy_constraint.append(RNA.energy_of_struct(current_seq, neg_structures[n]))
            pos_energy_constraint = []
            for p in range(0, len(pos_constraint)):
                pos_energy_constraint.append(RNA.energy_of_struct(current_seq, pos_constraint[p]))
            
            for x in range(0, len(pos_energy_constraint)):
                for k in range(0, len(neg_energy_constraint)):
                    if float(neg_energy_constraint[k]) - float(pos_energy_constraint[x]) < 0:
                        perfect = False
                        break

                if perfect:
                    break

            if perfect:
                break
        i += 1

    return Result(result_seq, score, pos_constraint, i) # result_seq

# objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)
def calculate_objective(sequence, structures):
    eos = []
    for struct in structures:
        eos.append(RNA.energy_of_struct(sequence, struct))
    
    gibbs = RNA.pf_fold(sequence)
    
    objective_difference_part = 0
    for i, value in enumerate(eos):
        for j in eos[i+1:]:
            objective_difference_part += math.fabs(value - j)
    
    return sum(eos) - len(eos) * gibbs[1] + 1 * objective_difference_part
    current_score = 0
    (current_seq_struct, current_seq_mfe) = RNA.fold(current_seq) # wird jetzt zwei mal berechnet, da oben nur in schleife, wenn struct in pos_constraints

    for l in range(0, len(pos_constraint)):
        eos_pos_cons = RNA.energy_of_struct(current_seq, pos_constraint[l])
        diff_eos_mfe = eos_pos_cons - current_seq_mfe
        current_score += diff_eos_mfe

    return current_score


if __name__ == "__main__":
    main()
