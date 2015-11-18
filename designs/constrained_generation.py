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
    parser.add_argument("-o", "--optimization", type=int, default=25000, help='Number of optimization iterations') # nicht mehr noetig
    parser.add_argument("-c", "--constraint", type=int, default=False, help='Number of negative constraints to record') # TODO: Name!
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

    if args.constraint:
        neg_structures = collections.deque(maxlen=args.constraint)
    else:
        if len(pos_constraint) < 3:
            neg_structures = collections.deque(maxlen=100)
        else:
            neg_structures = collections.deque(maxlen=20)

    # main loop from zero to number of solutions
    for i in range(0, args.number):
        r = optimization(dg, pos_constraint, neg_structures, args)
        if r is not None:
            r.write_out()
        else:
            sys.stdout.write("\r                                                             \r")
            sys.stdout.flush()
            print("No sequence with positive constraint.")


def optimization(dg, pos_constraint, neg_structures, args):
    score = float('Infinity')
    count = 0

    # randomly sample a initial sequence
    dg.sample()
    current_seq = dg.get_sequence()

    # number of mutations
    i = 0

    #while True:  # TODO progress passt so nicht mehr
    for i in range(0, args.optimization):

        (current_seq_struct, current_seq_mfe) = RNA.fold(current_seq)

        if current_seq_struct in pos_constraint:

            current_score = calculate_objective(current_seq, pos_constraint)

            if current_score < score:
                score = current_score
                result_seq = current_seq
                count = 0

                #if score == 0:
                   # return Result(result_seq, score, pos_constraint)

            else:
                count += 1
                if count > args.exit:
                    # return Result(result_seq, score, pos_constraint)
                    break

        if current_seq_struct not in neg_structures:
            neg_structures.append(current_seq_struct)

        while True:
            mut_nos = dg.sample_global()
            current_seq = dg.get_sequence()
            perfect = True

            if args.progress:
                sys.stdout.write("\rMutate global: {0:5.0f}/{1:5.0f} from NOS:{2:7.0f}".format(i, count, mut_nos))
                sys.stdout.flush()

            for x in range(0, len(pos_constraint)):
                pos_energy_constraint = RNA.energy_of_struct(current_seq, pos_constraint[x])

                for k in range(0, len(neg_structures)):

                    neg_energy_constraint = RNA.energy_of_struct(current_seq, neg_structures[k])

                    if float(neg_energy_constraint) - float(pos_energy_constraint) < 0:
                        perfect = False
                        break

                if perfect:
                    break

            if perfect:
                break
        i += 1

    return Result(result_seq, score, pos_constraint) # result_seq


def calculate_objective(current_seq, pos_constraint):
    current_score = 0
    (current_seq_struct, current_seq_mfe) = RNA.fold(current_seq) # wird jetzt zwei mal berechnet, da oben nur in schleife, wenn struct in pos_constraints

    for l in range(0, len(pos_constraint)):
        eos_pos_cons = RNA.energy_of_struct(current_seq, pos_constraint[l])
        diff_eos_mfe = eos_pos_cons - current_seq_mfe
        current_score += diff_eos_mfe

    return current_score


if __name__ == "__main__":
    main()
