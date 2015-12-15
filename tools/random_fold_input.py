from __future__ import print_function
import RNAdesign as rd
import RNA
import argparse


def get_structures(args):
    # calculate a RNAfold for random structures
    dg = rd.DependencyGraphMT(["." * args.length])
    structures = []
    amount = args.structures
    while (amount > 0):
        dg.sample()
        structure = RNA.fold(dg.get_sequence())[0]
        if (structure == "." * args.length):
            continue
        structures.append(structure)
        bipartite = rd.graph_is_bipartite(structures)
        if (bipartite):
            amount -= 1
        else:
            structures.pop()
    if (args.debug):
        print(*structures, sep="\n")
    # return structurs
    return structures

def main():
    parser = argparse.ArgumentParser(description='Generate a random sequence, get suboptimal structures and see if dependencygraph exists.')
    parser.add_argument("-f", "--filename", type=str, default=None, help='File name part of inp file to generate')
    parser.add_argument("-n", "--number", type=int, default=100, help='Number of input files to generate')
    parser.add_argument("-l", "--length", type=int, default=100, help='Length of random sequence')
    parser.add_argument("-s", "--structures", type=int, default=4, help='Number of structures as constraints input')
    parser.add_argument("-t", "--temperature", type=int, default=0, help='Temperature for the random sequence MFE folds (lower means more basepairs)')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()
    
    rd.initialize_library(args.debug)
    RNA.temperature = args.temperature
    
    for n in range(0, args.number):        
        # get random subopt structures
        structures = get_structures(args)
        if (args.filename is not None):
            fname = ".".join([args.filename, str(n), "inp"])
            with open(fname, 'w+') as f:
                f.write("\n".join(structures))
                f.write("\n\n;")
        else:
            print(">" + str(n))
            print("\n".join(structures))
        
if __name__ == "__main__":
    main()

