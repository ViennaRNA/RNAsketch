# correct_subopt.py
from __future__ import print_function

import re
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description='Takes the old RNAsubopt output and a motifconstraint to correct the energies.')
    parser.add_argument("-f", "--file", type = str, help='Specify RNAsubopt output.')
    parser.add_argument("-m", "--motif", type = str, help='Specify motif as string, e.g. --motif="GAUACCAG,(...((((,CCCUUGGCAGC,)...)))...),-9.22". The string is split at "," and [0] = sequence motif 5prime, [1] structure motif 5prime, [2] sequence motif 3prime, structure motif 3prime and [4] the bonus energy. Default is the example string which corresponds to the theophylline binding side.')
    parser.add_argument("-c", "--corOnly", default=False, action='store_true', help='Print only lines that have been corrected')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    args = parser.parse_args()

    sequence = ''
    structure = ''
    motif = args.motif.split(",") if (args.motif) else ["GAUACCAG", "(...((((", "CCCUUGGCAGC", ")...)))...)", -9.22]
    motif[4] = float(motif[4])
    f = open(args.file, "r") if (args.file)  else sys.stdin
    for line in f:
        line = line.rstrip('\r\n')
        a = line.split()
        if re.match(re.compile("^[\(\)\.\{\}\[\]\<\>\+\&]+$"), a[0], flags=0):
            if(sequence is ''):
                print("Could not find any sequence!", file=sys.stderr)
                exit
            else:
                structure = a[0]
                index1 = index2 = indexS1 = indexS2 = -1
                try:
                    index1 = sequence.index(motif[0])
                    index2 = sequence.index(motif[2])
                except ValueError, e:
                    pass
                if(index1 != -1 and index2 != -1):
                    try:
                        indexS1 = structure.index(motif[1])
                        indexS2 = structure.index(motif[3])
                    except ValueError, e:
                        pass
                    if(indexS1 == index1 and indexS2 == index2):
                        print(a[0] + "\t" + str(float(a[1])+motif[4]))
                    else:
                        if(not args.corOnly):
                            print(line)
        elif re.match(re.compile("^[ACGTU\&\+]+$"), a[0], flags=0) and sequence == '':
            sequence = a[0]
            print(line)
    f.close()
    

if __name__ == "__main__":
    main()

# End of file
