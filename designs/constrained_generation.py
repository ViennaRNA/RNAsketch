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

        print("Result...")
        print(self.sequence)
        print(self.mfe_struct)
        print(self.mfe_energy)

    '''def write_out(self):

        print(self.sequence)
        print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))'''


def main():
    parser = argparse.ArgumentParser(description='constrained sequence generation')
    parser.add_argument("-o", "--optimization", type=int, default=25000, help='Number of optimization iterations')
    parser.add_argument("-e", "--early_exit", type=int, default=10000, help='Exit optimization run if no better solution is aquired after early-exit trials.')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-i", "--input", default=False, action='store_true',
                        help='Read custom structures and sequence constraints from stdin')
    args = parser.parse_args()

    seq_constraint = ''
    pos_constraint = []
    # define structures
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
        seq_constraint = 'NNNNNNNNNNNNNNNNNNNNNNNNN'

    print("Positive constraint: "),
    print pos_constraint[0]

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

    check_constraint(dg, seq, computed_struct, pos_constraint, neg_structures, args)


def check_constraint(dg, sequence, structure, pos_constraint, neg_structures, args):
    if structure[0] != pos_constraint[0]:
        optimization(dg, structure, pos_constraint, neg_structures, args)

    else:
        return Result(sequence, structure)


def optimization(dg, structure, pos_constraint, neg_structures, args):
    count = 0
    result_seq = dg.get_sequence()
    result_struct = structure

    result_energy = RNA.energy_of_struct(result_seq, result_struct[0])

    for i in range(0, args.optimization):
        mut_nos = dg.mutate_global()
        new_seq = dg.get_sequence()
        (new_seq_struct, new_seq_mfe) = RNA.fold(new_seq)
        new_seq_energy = RNA.energy_of_struct(new_seq, new_seq_struct)

        if args.progress:
            sys.stdout.write("\rMutate global: {0:7.0f}/{1:5.0f} from NOS: {2:7.0f}".format(i, count, mut_nos))
            sys.stdout.flush()

        if float(result_energy) <= float(new_seq_energy) and new_seq_struct not in neg_structures:
            count += 1

            neg_structures.append(new_seq_struct)

            if count > args.early_exit:
                if new_seq_struct not in neg_structures:  # and float(new_seq_energy) < float(neg_structures_energy[-1])
                    neg_structures.append(new_seq_struct)
                break

        else:
            count = 0
            result_seq = new_seq
            result_energy = new_seq_energy
            result_struct = new_seq_struct

            if new_seq_struct not in neg_structures:
                if new_seq_struct == pos_constraint[0]:
                    print 'found'
                    break
                else:
                    neg_structures.append(new_seq_struct)

    return Result(result_seq, result_struct)


if __name__ == "__main__":
    main()
