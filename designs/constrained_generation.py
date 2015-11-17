import RNAdesign as rd
import RNA
import argparse
import sys
import collections
from collections import deque
import re


class Result:
	def __init__(self, sequence, score, pos_constraint):
		self.sequence = sequence
		self.pos_constraint = pos_constraint
		self.score = score
		self.eos = []

		(self.mfe_struct, self.mfe_energy) = RNA.fold(self.sequence)
		for struct in self.pos_constraint:
			self.eos.append(RNA.energy_of_struct(self.sequence, struct))		

	def write_out(self):
		# first clean up last line
		sys.stdout.write("\r                                                             \r")
		sys.stdout.flush()
		print(self.sequence + '\t{0:9.4f}'.format(self.score))
		for i, struct in enumerate(self.pos_constraint):
			print(struct + '\t{0:9.4f}\t{1:+9.4f}'.format(self.eos[i], self.eos[i]-self.mfe_energy))
		print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))


def main():
    parser = argparse.ArgumentParser(description='constrained sequence generation')
    parser.add_argument("-n", "--number", type=int, default=1, help='Number of designs to generate')
    parser.add_argument("-o", "--optimization", type=int, default=25000, help='Number of optimization iterations')
    parser.add_argument("-e", "--exit", type=int, default=1000, help='Exit optimization run if no better solution is aquired after exit trials.')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-i", "--input", default=False, action='store_true',
                        help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()
    rd.initialize_library(args.debug)
    
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
        # pos_constraint = ['((((....))))...(((....)))']
        pos_constraint = ['((((....))))....((((....))))........',
            '........((((....((((....))))....))))','((((((((....))))((((....))))....))))']
        seq_constraint = ''

	print("Positive constraint: ")
	print("\n".join(pos_constraint) + "\n" + seq_constraint)

    # construct dependency graph with these structures
    try:
        dg = rd.DependencyGraphMT(pos_constraint, seq_constraint)
    except Exception as e:
        print e
        quit()
	
    # collect sequences in order to exclude duplicates
    sequences = []

    neg_structures = collections.deque(maxlen=200/len(pos_constraint))
	
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

	score = float('Infinity')
	count = 0
    # randomly sample a initial sequence
	dg.set_sequence()
	current_seq = dg.get_sequence()

     # number of mutations
	for i in range(0, args.optimization):
		if args.progress:
			sys.stdout.write("\rOptimization iteration: {0:5.0f}/{1:5.0f}".format(i, args.optimization))
			sys.stdout.flush()

		
		(current_seq_struct, current_seq_mfe) = RNA.fold(current_seq)
     	
		if current_seq_struct in pos_constraint:
		    # calculate difference: sum all pos const (eos(pos) - mfe energy) -> should be smaller than previous solution
		    # if 0 or not better after args.exit trials return result, otherwise remember solution and go on.
			current_score = 0		

			for l in range(0, len(pos_constraint)):
				eos_pos_cons = RNA.energy_of_struct(current_seq, pos_constraint[l])
				eos_mfe = eos_pos_cons - current_seq_mfe 
				current_score += eos_mfe

			if current_score < score:
				score = current_score
				result_seq = current_seq
				count = 0
				if score == 0:
					return Result(result_seq, score, pos_constraint)
				
			else:
				count += 1
				if count > args.exit: 
					return Result(result_seq, score, pos_constraint) 
			
			# 
		
		if current_seq_struct not in neg_structures:
			neg_structures.append(current_seq_struct)

		#print neg_structures
		# neg_structures.append(current_seq_struct)
		while (True):
			dg.mutate_global()
			current_seq = dg.get_sequence()
			perfect = True
	
			for x in range(0, len(pos_constraint)):
				pos_energy_constraint = RNA.energy_of_struct(current_seq, pos_constraint[x])
				
				# print current_seq
				for k in range(0, len(neg_structures)):
					
					neg_energy_constraint = RNA.energy_of_struct(current_seq, neg_structures[k])
				
					if (float(neg_energy_constraint) - float(pos_energy_constraint) < 0):
						perfect = False
						break	
				
				if perfect:
					break
			
			if perfect:
				break		

	return Result(result_seq, score, pos_constraint) # result_seq


if __name__ == "__main__":
    main()
