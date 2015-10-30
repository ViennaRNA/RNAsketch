import RNAdesign as rd
import RNA
import argparse
import sys
import collections
from collections import deque

"""
class Result:
    def __init__(self, sequence, score, structures):
        self.sequence = sequence
        self.score = score
        self.structures = structures
        
        (self.mfe_struct, self.mfe_energy) = RNA.fold(self.sequence)
       
    def write_out(self):
        #first clean up last line
        sys.stdout.write("\r                                                             \r")
        sys.stdout.flush()
        print(self.sequence)
        print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))


def main():
 parser = argparse.ArgumentParser(description='constrained sequence generation')
    parser.add_argument("-s", "--sequence", type=int, default=4, help='Number of designs to generate')
"""   

	

# define structures
structures = ['((((....)))).............']
constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNN'
# construct dependency graph with these structures
seed = int(2)
dg = rd.DependencyGraphMT(structures, constraint, seed)

# randomly sample a initial sequence
dg.set_sequence()

seq = dg.get_sequence()
# print this sequence
#print seq


(struct, mfe) = RNA.fold(seq)

print "initial struct"
print struct

print "initial mfe"
print mfe

pos_constraint = '((((....))))...(((....)))'

neg_structures = collections.deque(maxlen=10)
#neg_structures = []
count = 0

count_better = 0

def optimization(seq, struct):
	dg.mutate_global()
    	new_seq = dg.get_sequence()
	(new_seq_struct, new_seq_mfe) = RNA.fold(seq)
	
	new_seq_energy = RNA.energy_of_struct(new_seq, new_seq_struct)
	seq_energy = RNA.energy_of_struct(seq, struct)
	
	# compare energy of struct
	if seq_energy <= new_seq_energy:
		optimization(seq, struct) #von urspruenglicher Sequenz ausgehend?
			#optimization(new_seq, new_seq_struct)				

			# check if already in list
		if new_seq_struct not in neg_structures:
			neg_structures.append(new_seq_struct)	
	
		"""
		global count
		count += 1
		
		# temporaer 
		if count > 500:
			sys.exit()
		else:
			optimization(seq, struct) #von urspruenglicher Sequenz ausgehend?
			#optimization(new_seq, new_seq_struct)				

			# check if already in list
			if new_seq_struct not in neg_structures:
				neg_structures.append(new_seq_struct)"""
		
		
		
	else:
		global count_better
		count_better += 1
		
		print "better sequence"
		print seq
		print struct
		print mfe
		check_constraint(new_seq, new_seq_struct)
		if count_bett > 2:
			sys.exit()

def check_constraint(sequence, structure):
	if pos_constraint != struct:
		neg_structures.append(seq)
		#print neg_structures

		#optimieren, mutieren...
		optimization(seq, struct)

	else:
		print "found"
		print seq
		print struct

check_constraint(seq, struct)	




		
	

	
