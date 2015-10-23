import RNAdesign as rd
import RNA
import argparse
import sys

# a tri-stable example target. (optional comment)
# ((((....))))....((((....))))........
# ........((((....((((....))))....))))
# ((((((((....))))((((....))))....))))
# below follows a simple (and optional) sequence constraint.
# CKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNB
# objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)


def main():
    parser = argparse.ArgumentParser(description='Design a tri-stable example same to Hoehner 2013 paper.')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    args = parser.parse_args()

    # define structures
    structures = ['((((....))))....((((....))))........',
                  '........((((....((((....))))....))))',
                  '((((((((....))))((((....))))....))))']
    
    # construct dependency graph with these structures
    dg = rd.DependencyGraphMT(structures)
    
    print "\n".join(structures) + "\n"
    # print the amount of solutions
    print 'Maximal number of solutions: ' + str(dg.number_of_sequences())
    # print the amount of connected components
    number_of_components = dg.number_of_connected_components()
    print 'Number of Connected Components: ' + str(number_of_components)

    for i in range(0, number_of_components):
        print '[' + str(i) + ']' + str(dg.component_vertices(i))
    print ''

    # resulting sequence
    result_seq = ''
    
    for n in range(0, args.number):
        score = 0
        count = 0
        
        # randomly sample a initial sequence
        dg.set_sequence()
        # print this sequence with score
        score = calculate_objective(dg.get_sequence(), structures);
        #print dg.get_sequence() + '\t' + str(score)
        
        # mutate globally for 1000 times and print
        for i in range(0, 500000):
            # write progress
            sys.stdout.write("\rMutate global: {0:7.0f}/{1:5.0f} from NOS: {2:7.0f}".format(i, count, dg.mutate_global()))
            sys.stdout.flush()
            
            this_score = calculate_objective(dg.get_sequence(), structures);
            #print dg.get_sequence() + '\t' + str(this_score)
            
            if (this_score < score):
                score = this_score
                result_seq = dg.get_sequence()
                count = 0
            else:
                count += 1
                if count > 25000:
                    break
        
        # finally print the result
        print_result(result_seq, structures, score)


# objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)
def calculate_objective(sequence, structures):
    eos = [0]
    for struct in structures:
        eos.append(RNA.energy_of_struct(sequence, struct))
    gibbs = RNA.pf_fold(sequence)
    
    return eos[1]+eos[2]+eos[3] - 3 * gibbs[1] + 1 * (pow(eos[1]-eos[2], 2) + pow(eos[1]-eos[3], 2) + pow(eos[2]-eos[3], 2))

# function which prints the result nicely to the screen
def print_result(sequence, structures, score):
    sys.stdout.write("\r                                                       \r")
    sys.stdout.flush()
    print sequence + '\t{0:4.5f}'.format(score)
    (mfe_struct, mfe_energy) = RNA.fold(sequence)
    
    for struct in structures:
        eos = RNA.energy_of_struct(sequence, struct)
        print struct + '\t{0:4.5f}\t{1:4.5f}'.format(eos, mfe_energy-eos)
    
    print mfe_struct + '\t{0:4.5f}'.format(mfe_energy)


if __name__ == "__main__":
    main()


