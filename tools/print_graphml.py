# http://hadim.fr/pygraphml/usage.html
import RNAdesign
import argparse
import tempfile
import sys
import re
from pygraphml import Graph
from pygraphml import GraphMLParser

def main():
    parser = argparse.ArgumentParser(description='Display a graphml file with python.')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    args = parser.parse_args()
    
    # define filename
    filename = ""
    parser = GraphMLParser()

    if (args.input and (args.graphml is None)):
        structures = []
        constraint = ""
        for line in sys.stdin:
            if re.match(re.compile("[\(\)\.]"), line, flags=0):
                structures.append(line.rstrip('\n'))
            elif re.match(re.compile("[ACGTUWSMKRYBDHVN]"), line, flags=0):
                constraint = line.rstrip('\n')
            elif re.search(re.compile("@"), line, flags=0):
                break;
        # get graphml from library
        graphml = RNAdesign.structures_to_graphml(structures, constraint)
        
        # write out a temporary graphml file
        with tempfile.NamedTemporaryFile() as f:
            f.write(graphml + "\n")
            f.flush()
            g = parser.parse(f.name)
            f.close()
    else:
        g = parser.parse(args.graphml)
    
    # now show the graph as plot
    g.show()
    
if __name__ == "__main__":
    main()

