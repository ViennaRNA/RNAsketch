# http://hadim.fr/pygraphml/usage.html
import RNAblueprint
import argparse
import tempfile
import sys
import re
import igraph

def graph_from_input(structures, constraint):
    g = None
    # get graphml from library
    try:
        graphml = RNAblueprint.structures_to_graphml(structures, constraint)
    except Exception as e:
        print e
        quit()
    
    # write out a temporary graphml file
    with tempfile.NamedTemporaryFile() as f:
        f.write(graphml + "\n")
        f.flush()
        g = igraph.Graph.Read_GraphML(f=f.name)
        f.close()
    return g

def main():
    parser = argparse.ArgumentParser(description='Display a Dependency Graph with python given structural constraints or a graphml file.')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Read graphml file with the given filename.')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-o", "--output", type=str, default=None, help='Write graph to PNG file with the given filename.')
    args = parser.parse_args()
    
    # define filename
    filename = ""
    g = igraph.Graph()
    
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
        g = graph_from_input(structures, constraint)
    elif (args.graphml is not None):
        g = igraph.Graph.Read_GraphML(f=args.graphml)
    elif (args.file is not None):
        structures = []
        constraint = ""
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
        g = graph_from_input(structures, constraint)
    else:
        print("ERROR: Please specify the input using command line arguments!")
        quit()
    
    # now show the graph as plot
    g.vs["label"] = [int(i) for i in g.vs["name"]]
    layout = g.layout_fruchterman_reingold()
    bool_dict = ["white", "red"]
    ear_dict=["black"]
    max_ear=int(max(g.es["ear"]))
    palette=igraph.RainbowPalette(n=max_ear+1)
    for c in range(0, max_ear+1):
        ear_dict.append(palette.get(c))
    igraph.plot(g, args.output, layout = layout, vertex_size=30, edge_width=5, margin=20, bbox=(800, 800), edge_color=[ear_dict[int(ear)] for ear in g.es["ear"]], vertex_color=[bool_dict[special] for special in g.vs["special"]])
    
    
if __name__ == "__main__":
    main()

