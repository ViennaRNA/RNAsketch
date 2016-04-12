# correct_subopt.py
from __future__ import print_function

#try:
#    from PyDesign import *
#except ImportError, e:
#    print(e.message)
#    exit(1)
import PyDesign as pd
import re
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Takes the old RNAsubopt output and a softconstraint to correct the energies.')
    parser.add_argument("-f", "--file", type = str, help='Specify RNAsubopt output.')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    f = open(args.file, 'r')
    sequence = ''
    structure = ''
    ligand = ["GAUACCAG", "(...((((", "CCCUUGGCAGC", ")...)))...)", -9.22]
    for line in f:
        a = line.split()
        if re.match(re.compile("^[\(\)\.\{\}\[\]\<\>\+\&]+$"), a[0], flags=0):
            if(sequence is ''):
                print("Could not find any sequence!", file=sys.stderr)
                exit
            else:
                structure = a[0]
                index1 = index2 = indexS1 = indexS2 = -1
                try:
                    index1 = sequence.index(ligand[0])
                    index2 = sequence.index(ligand[2])
                except ValueError, e:
                    pass
                if(index1 != -1 and index2 != -1):
                    try:
                        indexS1 = structure.index(ligand[1])
                        indexS2 = structure.index(ligand[3])
                    except ValueError, e:
                        pass
                    if(indexS1 == index1 and indexS2 == index2):
                        print(a[0] + "\t" + str(float(a[1])+ligand[4]))
                    else:
                        print(line)
        elif re.match(re.compile("^[ACGTU\&\+]+$"), a[0], flags=0) and sequence == '':
            sequence = a[0]
    f.close()
   

if __name__ == "__main__":
    main()

# End of file
