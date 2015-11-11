import RNAdesign as rd
import RNA
import argparse
import sys
import collections
from collections import deque
import re


class Result:
    def __init__(self, sequence, structure, pos_constraint, sequences):
        self.sequence = sequence
        self.structure = structure
        self.pos_constraint = pos_constraint
        self.sequences = sequences

        (self.mfe_struct, self.mfe_energy) = RNA.fold(self.sequence)

    def write_out(self):
        # first clean up last line
        sys.stdout.write("\r                                                             \r")
        sys.stdout.flush()
        print(self.sequence + '\t{0:9.4f}'.format(self.mfe_energy))
        print(self.mfe_struct)


def main():
    parser = argparse.ArgumentParser(description='constrained sequence generation')
    parser.add_argument("-n", "--number", type=int, default=1, help='Number of designs to generate')
    parser.add_argument("-o", "--optimization", type=int, default=25000, help='Number of optimization iterations')
    parser.add_argument("-e", "--early_exit", type=int, default=10000, help='Exit optimization run if no better solution is aquired after early-exit trials.')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-i", "--input", default=False, action='store_true',
                        help='Read custom structures and sequence constraints from stdin')
    args = parser.parse_args()

    # define structures
    seq_constraint = ''
    pos_constraint = []
    if args.input:
        for line in sys.stdin:
            if re.match(re.compile("[\(\)\.]"), line, flags=0):
                pos_constraint.append(line.rstrip('\n'))
            elif re.match(re.compile("[ACGTUWSMKRYBDHVN]"), line, flags=0):
                seq_constraint = line.rstrip('\n')
            elif re.search(re.compile("@"), line, flags=0):
                break

    else:
        pos_constraint = ['((((....))))...(((....)))']
        seq_constraint = ''

    print("Positive constraint: ")
    print(pos_constraint[0])

    if seq_constraint is not '':
        print("Sequence constraint: ")
        print(seq_constraint)
    print

    # construct dependency graph with these structures
    try:
        dg = rd.DependencyGraphMT(pos_constraint, seq_constraint)
    except Exception as e:
        print e
        quit()
	
    # collect sequences in order to exclude duplicates
    sequences = []

    neg_structures = collections.deque(maxlen=100)
	
    # main loop from zero to number of solutions
    for i in range(0, args.number):
        r = optimization(dg, pos_constraint, neg_structures, sequences, args)
        if r is not None:
            r.write_out()
        else:
            sys.stdout.write("\r                                                             \r")
            sys.stdout.flush()
            print("No sequence with positive constraint.")


def optimization(dg, pos_constraint, neg_structures, sequences, args):
    # randomly sample a initial sequence
	dg.set_sequence()
	result_seq = dg.get_sequence()

     # number of mutations
	for i in range(0, args.optimization):
		if args.progress:
			sys.stdout.write("\rOptimization iteration: {0:5.0f}/{1:5.0f}".format(i, args.optimization))
			sys.stdout.flush()

		
		(result_seq_struct, result_seq_mfe) = RNA.fold(result_seq)
     	
		if result_seq_struct == pos_constraint[0]:
			return Result(result_seq, result_seq_struct, pos_constraint, sequences)
		
		if result_seq_struct not in neg_structures:
			neg_structures.append(result_seq_struct)

		print neg_structures
		# neg_structures.append(result_seq_struct)
		while (True):
			dg.mutate_global()
			result_seq = dg.get_sequence()
			pos_energy_constraint = RNA.energy_of_struct(result_seq, pos_constraint[0])
			perfect = True
			# print result_seq
			for k in range(0, len(neg_structures)):
				neg_energy_constraint = RNA.energy_of_struct(result_seq, neg_structures[k])
				
				if float(pos_energy_constraint) > float(neg_energy_constraint):
					# print pos_energy_constraint
					# print neg_energy_constraint
					perfect = False
					break
				
			if perfect:
				break
			
					

	return Result(result_seq, result_seq_struct, pos_constraint, sequences)


if __name__ == "__main__":
    main()
