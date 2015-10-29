import RNAdesign as rd
import RNA

# define structures
structures = ['((((....)))).............']
constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNN'
# construct dependency graph with these structures
dg = rd.DependencyGraphMT(structures, constraint)

# randomly sample a initial sequence
dg.set_sequence()

seq = dg.get_sequence()
# print this sequence
print seq



(struct, mfe) = RNA.fold(seq)
print struct
#print mfe

pos_constraint = '((((....))))...(((....)))'

neg_structures = []

if pos_constraint == struct:
	neg_structures.append(seq)
	
	#optimieren, mutieren...
	
else:
	print seq
	print struct

def optimization(seq):
	dg.mutate_local()
	




