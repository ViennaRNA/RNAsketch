import RNAdesign as rd
import RNA
import argparse
import sys
import collections
from collections import deque
import re


class Result:
    def __init__(self, sequence, structure):
        self.sequence = sequence
        self.structure = structure

        (self.mfe_struct, self.mfe_energy) = RNA.fold(self.sequence)

        print("temporary result...")
        print(self.sequence)
        print(self.mfe_struct)
        print(self.mfe_energy)

    '''def write_out(self):

        print(self.sequence)
        print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))'''


def main():

    '''parser = argparse.ArgumentParser(description='constrained sequence generation')
    parser.add_argument("--structure", help='Specify structure')
    args = parser.parse_args()
    print args.structure


    # define structures
    pos_constraint = []
    if (args.structure):

        pos_constraint = args.structure

    else:'''
    pos_constraint = ['((((....))))...(((....)))']
    seq_constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNN'

    print("Positive constraint: "),
    print pos_constraint[0]

    global count
    count = 0
    # construct dependency graph with these structures
    seed = int(2)
    dg = rd.DependencyGraphMT(pos_constraint, seq_constraint, seed)

    # randomly sample a initial sequence
    dg.set_sequence()

    seq = dg.get_sequence()

    (struct, mfe) = RNA.fold(seq)
    computed_struct = [struct]

    print "initial struct: ",
    print struct
    print "initial seq: ",
    print seq

    print "initial mfe",
    print mfe
    neg_structures = collections.deque(maxlen=10)
    neg_structures_energy = collections.deque(maxlen=10)

    check_constraint(dg, seq, computed_struct, pos_constraint, neg_structures, neg_structures_energy)


def check_constraint(dg, sequence, structure, pos_constraint, neg_structures, neg_structures_energy):
    #print 'checking constraints'

    if structure[0] != pos_constraint[0]:
        neg_energy = RNA.energy_of_struct(sequence, structure[0])

        neg_structures.append(structure[0])
        #neg_structures_energy.append(neg_energy)
        # optimieren, mutieren... #neuen dg?
        optimization(dg, sequence, structure, pos_constraint, neg_structures, neg_structures_energy)

    else:
        return Result(sequence, structure)


def optimization(dg, seq, struct, pos_constraint, neg_structures, neg_structures_energy):
    #print "optimizing"
    global count
    dg.mutate_global()

    old_seq = seq
    old_struct = struct

    new_seq = dg.get_sequence()
    (new_seq_struct, new_seq_mfe) = RNA.fold(new_seq)

    new_seq_struct_list = [new_seq_struct]

    new_seq_energy = RNA.energy_of_struct(new_seq, new_seq_struct)
    seq_energy = RNA.energy_of_struct(old_seq, struct[0])

    if float(seq_energy) <= float(new_seq_energy):
        count += 1
        if count > 100:
            print "Die letzten hundert Mutationen fuehren zu keinem besseren Ergebnis."
            print neg_structures
            sys.exit(Result(old_seq, old_struct))

        # optimize nur, wenn structure nicht schon in negative constraint list
        if new_seq_struct not in neg_structures: # and float(new_seq_energy) < float(neg_structures_energy[-1]):
            neg_structures.append(new_seq_struct)

        optimization(dg, old_seq, struct, pos_constraint, neg_structures, neg_structures_energy)

    else:
        count = 0
        dg = rd.DependencyGraphMT(pos_constraint, new_seq)
        neg_structures.append(old_struct[0])
        check_constraint(dg, new_seq, new_seq_struct_list, pos_constraint, neg_structures, neg_structures_energy)

if __name__ == "__main__":
    main()

