from __future__ import print_function

try:
    from PyDesign import *
except ImportError, e:
    print(e.message)
    exit(1)

import RNAdesign as rd
import argparse
import sys
import time
from pylab import remainder

def main():
    parser = argparse.ArgumentParser(description='Design of a transcription regulating riboswitch similar to Wachsmuth et al. 2013.')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-q", "--nupack", default=False, action='store_true', help='Use Nupack instead of the ViennaRNA package (for pseudoknots)')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    parser.add_argument("-e", "--exit", type=int, default=500, help='Exit optimization run if no better solution is aquired after (exit) trials.')
    parser.add_argument("-m", "--mode", type=str, default='random', help='Mode for getting a new sequence: sample, sample_local, sample_global, random')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    #### new
    parser.add_argument("-r", "--range", type=str, default='6,20', help='Range of which a random number is generated to determine a spacer length')
    parser.add_argument("-R", "--ratio", type=str, default='9:1', help='Ratio of the alternative and binding competent state')
    parser.add_argument("-s", "--start", type=int, default=26, help='Position within the aptamer from which on the alternative state should be made')
    parser.add_argument("-l", "--altl", type=int, default=11, help='Length of the stem in the alternative state')
    
    args = parser.parse_args()

    print("# Options: number={0:d}, exit={1:d}, mode={2:}, nupack={3:}".format(args.number, args.exit,  args.mode, str(args.nupack)))
    rd.initialize_library(args.debug, args.kill)

    try:
        spacerrange = map(int, args.range.split(",")) #convert region to array of integers
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)
        sys.exit()

    try:
        ratio = map(int, args.ratio.split(":")) #convert ratio to array of integers
    except Exception as e:
        print( "Error: %s" % e , file=sys.stderr)
        sys.exit()
        
    if not (sum(ratio) <= 10):
        sys.stderr.write("You specified ratios that in total become larger than 1\n")
        sys.exit()

    for i, a in enumerate(ratio):
        ratio[i]=a/float(10)
        print(ratio[i])
        
    aptseq = ''
    aptstr = []
    if (args.input):
        data = ''
        for line in sys.stdin:
            data = data + '\n' + line
        (aptstr, aptseq) = read_input(data)
    elif (args.file is not None):
        print("# Input File: {0:}".format(args.file))
        (aptstr, aptseq) = read_inp_file(args.file)
    else:
        aptstr = ['(((((...((((((((.....)))))...)))...)))))..']
        aptseq = 'AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCA'

    if (((args.start-1)+args.altl) > len(aptseq)):
        sys.stderr.write("You specified -s "+str(args.start)+" and -l "+str(args.altl)+" that gets out of range of the aptamer length " + str(len(aptseq)) + "\n")
        sys.exit()
    aptposS = args.start
    aptposL = args.altl
    comp = int(len(aptstr[0])/2) # length of the sequence that could get complementary

    # main loop from zero to number of solutions
    for n in range(0, args.number):
        global_score = float("inf")
        global_design = None
       
        for spacer in range(spacerrange[0], (spacerrange[1]+1)):
            
            # Build input structures and sequence constraints
            structures = []
            constraint = ''
            start_sequence = ''
            print("# spacer length: " + str(spacer) + "\tlength of the complementary region: " + str(comp))

            ONL = aptstr[0] + ("."*spacer) + ("."*comp)
            ONH = aptstr[0] + ("."*spacer) + ("."*comp)
            
            OFFL = ("."*(aptposS-1)) + ("("*aptposL) + ("."*(len(aptseq)-(aptposS+aptposL))) + ("."*spacer) + (")"*aptposL) + ("."*(comp-aptposL+1))
            OFFH = ("."*(aptposS-1)) + ("<"*(len(aptseq)-aptposS+1)) + ("."*spacer) + ("."*comp)

            OFFI = ("."*len(aptstr[0])) + ("."*spacer) + ("."*comp) # simple state or come up with a usefull intermediat like "........((((((((.....)))))...)))((((...............)))).................."
            
            SEQ = aptseq + ("N"*(spacer+comp))
            print("# Used input:\n# " + ONL + "\n# " + OFFI + "\n# " + OFFL + "\n#\n# " + ONH + "\n# " + OFFI + "\n# " + OFFH + "\n# " + SEQ)



            # try to construct dependency graph, catch errors and timeouts
            dg = None
            construction_time = 0.0
            sample_time = 0.0
            
            # construct dependency graph with these structures
            try:
                start = time.clock()
                dg = rd.DependencyGraphMT([ONL, OFFI, OFFL], SEQ)
                construction_time = time.clock() - start
            except Exception as e:
                print( "Error: %s" % e , file=sys.stderr)
                
            # general DG values
            print("# " + "\n# ".join(structures) + "\n# " + SEQ)

            if (dg is not None):
                
                # if requested write out a graphml file
                if args.graphml is not None:
                    with open(args.graphml, 'w') as f:
                        f.write(dg.get_graphml() + "\n")
                        
                # print the amount of solutions
                print('# Maximal number of solutions: ' + str(dg.number_of_sequences()))
                # print the amount of connected components
                number_of_components = dg.number_of_connected_components()
                print('# Number of Connected Components: ' + str(number_of_components))
                for i in range(0, number_of_components):
                    print('# [' + str(i) + ']' + str(dg.component_vertices(i)))
                        
                # remember general DG values
                graph_properties = get_graph_properties(dg)

                if (args.nupack):
                    design = nupackDesign([], start_sequence)
                else:
                    design = vrnaDesign([], start_sequence)

                #[0]: Ligand bound state including hard constraint   = on + ligand + hardconstr1
                #[1]: OFF state with hard constraints of alternative = off + hardconstr2
                #[2]: OFF state with hard constraint of intermediate = off + hardconstr3
                #[3]: Ligand bound state without hard constraints    = on + ligand
                #[4]: OFF state without hard constraints             = off
                #[5]: OFF state with hard constraint of on           = on in off

                
                design.newState("on ligand h1", ONL)
                design.state["on ligand h1"].ligand=["GAUACCAG&CCCUUGGCAGC", "(...((((&)...)))...)", -9.22]
                design.state["on ligand h1"].constraint=ONH #design.state["on ligand h1"].enforce_constraint=True

                
                design.newState("off h2", OFFL)
                design.state["off h2"].constraint=OFFH

                design.newState("off h3", OFFI)
                design.state["off h3"].constraint=OFFI

                design.newState("on ligand", ONL)
                design.state["on ligand"].ligand=["GAUACCAG&CCCUUGGCAGC", "(...((((&)...)))...)", -9.22]

                design.newState("off", OFFL) # no structure just to calculate partition function of the sequence

                design.newState("on in off",OFFI)
                design.state["on in off"].constraint=ONH
                                
                    
                # print header for csv file
                if (args.csv):
                    print(";".join(["exit",
                                    "mode",
                                    "score",
                                    "num_mutations",
                                    "construction_time",
                                    "sample_time",
                                    design.write_csv_header()] +
                                   graph_properties.keys()))
                    
                #add variable
                design.aptstr=aptstr
                design.aptsp_length=len(aptseq)+spacer
                design.ratio=ratio
                
                start = time.clock()
                # do a random sampling
                (score, number_of_mutations) = classic_optimization(dg, design, objective_function=calculate_switch_objective, exit=args.exit, mode='random', progress=args.progress)
                # sum up for a complete number of mutations
                number_of_mutations += number_of_mutations
                sample_time = time.clock() - start

                #print the details of the objective function
                design.sequence = "AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAGUUGCCGGGGGACCGGCGGAUG"
                calculate_switch_objective(design, printDetails=True)
                if (args.csv):
                    print(args.exit,
                          "\"" + args.mode + "\"",
                          score,
                          number_of_mutations,
                          construction_time,
                          sample_time,
                          design.write_csv(),
                          *graph_properties.values(), sep=";")
                else:
                    print(design.write_out(score))
                    
                if(score<global_score):
                    global_score = score
                    global_design = design
            else:
                print('# Construction time out reached!')
            
    print("\n\n" + "Best desgin:\n" + global_design.write_out(global_score))

def calculate_switch_objective(design, printDetails=False):
    '''Calculates the objective function given a design object containing
    the partition function energies of the following four states
    [0]: Ligand bound state including hard constraint   = on + ligand + hardconstr1
    [1]: OFF state with hard constraints of alternative = off + hardconstr2
    [2]: OFF state with hard constraint of intermediate = off + hardconstr3
    [3]: Ligand bound state without hard constraints    = on + ligand
    [4]: OFF state without hard constraints             = off
    [5]: OFF state with hard constraint of on           = on in off

    
    Objective function two state:   ([0] - [2]) + ([1] - [3])
    Objective function three state: ([0] - [2]) + ([1] - [3]) + ([4] - [3])

    Objective function prob state:  (1 - P(A|+)) + |(P(0|-)-0.9)| + |(P(A|-)-0.1)|
    Objective function prob state:  (1 - |[3]-[0]|) + (|[4]-[1]|-0.9) + (|[4]-[5]|-0.1)
    :param design: Design object containing all necessary partition function energies

    '''
    
    term1 = (1-Z_from_G((design.state["on ligand h1"].pf_energy-design.state["on ligand"].pf_energy))) #(A|+)
    if(printDetails): print("dG=%.3f Z=%.3f term1=%.3f" % ((design.state["on ligand h1"].pf_energy-design.state["on ligand"].pf_energy), Z_from_G((design.state["on ligand h1"].pf_energy-design.state["on ligand"].pf_energy)), term1))
    
    term2 = abs(Z_from_G((design.state["off h2"].pf_energy-design.state["off"].pf_energy))-design.ratio[0]) #(0|-)
    if(printDetails): print("dG=%.3f Z=%.3f R=%.3f term2=%.3f" % ((design.state["off h2"].pf_energy-design.state["off"].pf_energy), Z_from_G((design.state["off h2"].pf_energy-design.state["off"].pf_energy)), design.ratio[0], term2))
    
    term3 = abs(Z_from_G((design.state["on in off"].pf_energy-design.state["off"].pf_energy))-design.ratio[1]) #(A|-)
    if(printDetails): print("dG=%.3f Z=%.3f R=%.3f term3=%.3f" % ((design.state["on in off"].pf_energy-design.state["off"].pf_energy), Z_from_G((design.state["on in off"].pf_energy-design.state["off"].pf_energy)), design.ratio[1], term3))
    

    if(printDetails): print("(A|+)=%.3f (0|-)=%.3f (A|-)=%.3f" % (term1, term2, term3))

    return (term1 + term2 + term3)

def calculate_bp_distance(s1, s2):
    '''
    Compare two structures of unequal length and determine the base pair distance 
    d_BP(P_a,P_b) = |BP_a| + |BP_b| - 2|BP_a \cap  BP_b| Note that we take only the 
    sub string of the longer sequence equal to the length of the shorter one.
    :param s1: First structure to be compared
    :param s2: Second structure to be compared
    '''

    #re-order and make s1 the shorter string
    if(len(s1)>len(s2)):
        t=s1
        s1=s2
        s2=s1
        
    bpt1=create_bp_table(s1)
    bpt2=create_bp_table(s2)
    inboth=0
    for i, val in enumerate(bpt1):
        if(val!=-1 and bpt1[i]==bpt2[i]):
            inboth+=1
    dist = s1.count(")") + s2[0:len(s1)].count(")") - 2 * inboth
    return (inboth, dist)

def G_from_Z(Z):
    return - ((37.0 + 273.15)*1.98717)/1000.0 * math.log(Z)

def Z_from_G(G):
    return math.exp(- (G/ (((37.0 + 273.15)*1.98717)/1000.0)))


if __name__ == "__main__":
    main()




