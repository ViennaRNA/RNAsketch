# http://hadim.fr/pygraphml/usage.html
import argparse
from pygraphml import Graph
from pygraphml import GraphMLParser

def main():
    parser = argparse.ArgumentParser(description='Display a graphml file with python.')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    args = parser.parse_args()
    
    parser = GraphMLParser()
    g = parser.parse(args.graphml)
    g.show()
    
if __name__ == "__main__":
    main()

