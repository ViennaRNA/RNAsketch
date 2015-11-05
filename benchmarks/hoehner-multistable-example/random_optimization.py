from __future__ import print_function
import RNAdesign as rd
import RNA
import argparse
import sys
import re

# a tri-stable example target. (optional comment)
# ((((....))))....((((....))))........
# ........((((....((((....))))....))))
# ((((((((....))))((((....))))....))))
# below follows a simple (and optional) sequence constraint.
# CKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNB
# objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)

class Result:
    def __init__(self, sequence, score, structures):
        self.sequence = sequence
        self.score = score
        self.structures = structures
        self.eos = []
        
        (self.mfe_struct, self.mfe_energy) = RNA.fold(self.sequence)
        for struct in self.structures:
            self.eos.append(RNA.energy_of_struct(self.sequence, struct))
    def write_out(self):
        #first clean up last line
        sys.stdout.write("\r" + " " * 60 + "\r")
        sys.stdout.flush()
        print(self.sequence + '\t{0:9.4f}'.format(self.score))
        for i, struct in enumerate(self.structures):
            print(struct + '\t{0:9.4f}\t{1:+9.4f}'.format(self.eos[i], self.eos[i]-self.mfe_energy))
        print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))

def main():
    parser = argparse.ArgumentParser(description='Design a tri-stable example same to Hoehner 2013 paper.')
    parser.add_argument("-n", "--number", type=int, default=40, help='Number of designs to generate')
    parser.add_argument("-x", "--optimization_min", type=int, default=0, help='Minimal number of random optimization iterations')
    parser.add_argument("-y", "--optimization_max", type=int, default=50000, help='Maximal number of random optimization iterations')
    parser.add_argument("-m", "--mode", type=str, default='sample', help='Mode for getting a new sequence: sample, mutate_local, mutate_global')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    args = parser.parse_args()
    
    print ("# random_optimization.py")
    print ("# Options: number={0:d}, optimization_min={1:d}, optimization_max={2:d}, mode={3:}".format(args.number, args.optimization_min, args.optimization_max, args.mode))

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
    else:
        structures = ['((((....))))....((((....))))........',
            '........((((....((((....))))....))))',
            '((((((((....))))((((....))))....))))']
        constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
    
    # construct dependency graph with these structures
    dg = rd.DependencyGraphMT(structures, constraint)
    
    print("# " + "\n# ".join(structures) + "\n# " + constraint)
    # print the amount of solutions
    print('# Maximal number of solutions: ' + str(dg.number_of_sequences()))
    # print the amount of connected components
    number_of_components = dg.number_of_connected_components()
    print('# Number of Connected Components: ' + str(number_of_components))

    for i in range(0, number_of_components):
        print('# [' + str(i) + ']' + str(dg.component_vertices(i)))
    
    
    # optmizations start here
    for optimization_iterations in xrange(args.optimization_min, args.optimization_max, 500):
        for n in range(0, args.number):
            r = optimization_run(dg, structures, args, optimization_iterations)
            
            # process result and write result of this optimization to stdout
            eos_diff = []
            reached_mfe = 0
            for i in range(0, len(r.structures)):
                eos_diff.append(r.eos[i]-r.mfe_energy)
                if (r.structures[i] == r.mfe_struct):
                    reached_mfe = 1
            
            if (args.progress):
                sys.stdout.write("\r" + " " * 60 + "\r")
                sys.stdout.flush()
            print (optimization_iterations, r.score, reached_mfe, *eos_diff, sep=";")

# main optimization
def optimization_run(dg, structures, args, optimization_iterations):
    score = 0
    
    # randomly sample a initial sequence
    dg.set_sequence()
    result_sequence = dg.get_sequence()
    # print this sequence with score
    score = calculate_objective(dg.get_sequence(), structures);
    #print dg.get_sequence() + '\t' + str(score)
    
    # mutate globally for num_opt times and print
    for i in range(0, optimization_iterations):
        # mutate sequence completely random
        if args.mode == 'sample':
            mut_nos = dg.set_sequence()
        elif args.mode == 'mutate_global':
            mut_nos = dg.mutate_global()
        elif args.mode == 'mutate_local':
            mut_nos = dg.mutate_local()
        else:
            sys.stdout.write("Wrong mutate argument: " + args.mutate + "\n")
            sys.exit(1)
        # write progress
        if (args.progress):
            sys.stdout.write("\rMutate global: {0:7.0f} from NOS: {1:7.0f}".format(i, mut_nos) + " " * 20)
            sys.stdout.flush()
        
        this_score = calculate_objective(dg.get_sequence(), structures);
        
        if (this_score < score):
            score = this_score
            result_sequence = dg.get_sequence()
    
    # finally return the result
    return Result(result_sequence, score, structures)

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


