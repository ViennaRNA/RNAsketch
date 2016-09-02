# http://hadim.fr/pygraphml/usage.html
import RNAblueprint
import argparse
import tempfile
import sys
import re

def main():
    parser = argparse.ArgumentParser(description='Generate a graphml file with python given structural constraints.')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-o", "--output", type=str, default=None, help='Write graphml file with the given filename.')
    args = parser.parse_args()
    
    # define filename
    structures = []
    constraint = ""
    
    if (args.input):
        for line in sys.stdin:
            if re.match(re.compile("[\(\)\.]"), line, flags=0):
                structures.append(line.rstrip('\n'))
            elif re.match(re.compile("[ACGTUWSMKRYBDHVN]"), line, flags=0):
                constraint = line.rstrip('\n')
            elif re.search(re.compile("@"), line, flags=0):
                break;
    elif (args.file is not None):
        with open(args.file) as f:
            data = f.read()
            lines = data.split("\n")
            for line in lines:
                if re.match(re.compile("[\(\)\.]"), line):
                    structures.append(line)
                if re.match(re.compile("[\ AUGC]"), line):
                    elements = str(line)
                    constraint = elements.replace(" ", "N")
                if line.startswith(";"):
                    break
    else:
        print("ERROR: Please specify the input using command line arguments!")
        quit()
    
    try:
        graphml = RNAblueprint.structures_to_graphml(structures, constraint)
    except Exception as e:
        print e
        quit()
    
    # write out a graphml file
    if (args.output is not None):
        with open(args.output, 'w+') as f:
            f.write(graphml + "\n")
            f.flush()
    else:
        print(graphml + "\n")
    
if __name__ == "__main__":
    main()

