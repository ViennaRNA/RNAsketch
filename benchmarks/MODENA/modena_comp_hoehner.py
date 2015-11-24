from __future__ import print_function
import RNAdesign as rd
import RNA
import argparse
import sys
import re
import math
from collections import deque
import collections

# a tri-stable example target. (optional comment)
# ((((....))))....((((....))))........
# ........((((....((((....))))....))))
# ((((((((....))))((((....))))....))))
# below follows a simple (and optional) sequence constraint.
# CKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNB
# objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)

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
        '''
        #first clean up last line
        sys.stdout.write("\r" + " " * 60 + "\r")
        sys.stdout.flush()
        print(self.sequence + '\t{0:9.4f}'.format(self.score))
        for i, struct in enumerate(self.structures):
            print(struct + '\t{0:9.4f}\t{1:+9.4f}\t{2:9.4f}'.format(self.eos[i], self.eos[i]-self.mfe_energy, self.probs[i]))
        print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))'''


def main():
    parser = argparse.ArgumentParser(description='Design a tri-stable example same to Hoehner 2013 paper.')
    parser.add_argument("-n", "--number", type=int, default=100, help='Number of designs to generate')
    parser.add_argument("-j", "--jump", type=int, default=300, help='Do random jumps in the solution space for the first (jump) trials.')
    parser.add_argument("-e", "--exit", type=int, default=500, help='Exit optimization run if no better solution is aquired after (exit) trials.')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-f", "--file", type = str, default=False, help='Read file in *.inp format')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    print ("# Options: file={0:s}, number={1:d}, jump={2:d}, exit={3:d}".format(args.file, args.number, args.jump, args.exit))
    rd.initialize_library(args.debug)
    # define structures
    structures = []
    constraint = ""

    if (args.file):
        filename = args.file
        # filename = "alpha_operon_edif0.0.inp"
        with open(filename) as f:
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

    
    # construct dependency graph with these structures
    try:
        dg = rd.DependencyGraphMT(structures, constraint)
    except Exception as e:
        print (e)
        quit()
    
    print("# " + "\n# ".join(structures) + "\n# " + constraint)
    # print the amount of solutions
    print('# Maximal number of solutions: ' + str(dg.number_of_sequences()))
    # print the amount of connected components
    number_of_components = dg.number_of_connected_components()
    print('# Number of Connected Components: ' + str(number_of_components))

    for i in range(0, number_of_components):
        print('# [' + str(i) + ']' + str(dg.component_vertices(i)))
    print('')
    
    # if requested write out a graphml file
    if args.graphml is not None:
        with open(args.graphml, 'w') as f:
            f.write(dg.get_graphml() + "\n")
    
    # main loop from zero to number of solutions

    # optmizations start here
    mfe_reached_str = ""
    diff_eos_mfe_str = ""
    for s in range(0, len(structures)):
        mfe_reached_str = mfe_reached_str + "mfe_reached_" + str(s) +";"
        diff_eos_mfe_str = diff_eos_mfe_str + "diff_eos_mfe_" + str(s) + ";"
    print ("num_mutations;seq_length;sequence;" + mfe_reached_str + diff_eos_mfe_str)
    for n in range(0, args.number):
        r = optimization_run(dg, structures, args)
        r.write_out()

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
     
        print (r.number_of_mutations, len(r.sequence), "\"" + r.sequence + "\"", *(mfe_reached + diff_eos_mfe), sep=";")

# main optimization
def optimization_run(dg, structures, args):
    score = 0
    count = 0
    jumps = args.jump
    # randomly sample a initial sequence
    dg.sample()
    # print this sequence with score
    score = calculate_objective(dg.get_sequence(), structures);
    #print dg.get_sequence() + '\t' + str(score)
    
    # sample globally for num_opt times and print
    i = 0
    while 1:
        # sample sequence
        if jumps:
            mut_nos = dg.sample()
            jumps -= 1
            count = 0
        else:
            mut_nos = dg.sample_global()
        # write progress
        if (args.progress):
            sys.stdout.write("\rMutate global: {0:7.0f}/{1:5.0f} from NOS: {2:7.0f}".format(i, count, mut_nos) + " " * 20)
            sys.stdout.flush()
        
        this_score = calculate_objective(dg.get_sequence(), structures);
        
        if (this_score < score):
            score = this_score
            count = 0
        else:
            dg.revert_sequence();
            count += 1
            if count > args.exit:
                break
        i += 1
    
    # finally return the result
    return Result(dg.get_sequence(), score, structures, i)

# objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)
def calculate_objective(sequence, structures):
    eos = []
    for struct in structures:
        eos.append(RNA.energy_of_struct(sequence, struct))
    
    gibbs = RNA.pf_fold(sequence)
    
    objective_difference_part = 0
    for i, value in enumerate(eos):
        for j in eos[i+1:]:
            objective_difference_part += pow(value - j, 2)
    
    return sum(eos) - len(eos) * gibbs[1] + 1 * objective_difference_part

if __name__ == "__main__":
    main()



