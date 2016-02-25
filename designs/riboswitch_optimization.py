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

def main():
    parser = argparse.ArgumentParser(description='Design of a transcription regulating riboswitch similar to Wachsmuth et al. 2013.')
    parser.add_argument("-f", "--file", type = str, default=None, help='Read file in *.inp format')
    parser.add_argument("-i", "--input", default=False, action='store_true', help='Read custom structures and sequence constraints from stdin')
    parser.add_argument("-q", "--nupack", default=False, action='store_true', help='Use Nupack instead of the ViennaRNA package (for pseudoknots)')
    parser.add_argument("-n", "--number", type=int, default=4, help='Number of designs to generate')
    parser.add_argument("-j", "--jump", type=int, default=300, help='Do random jumps in the solution space for the first (jump) trials.')
    parser.add_argument("-e", "--exit", type=int, default=500, help='Exit optimization run if no better solution is aquired after (exit) trials.')
    parser.add_argument("-s", "--strelem", type=int, default=1800, help='Optimize structural elements and exit after (strelem) unsucessful trials.')
    parser.add_argument("-m", "--mode", type=str, default='sample_global', help='Mode for getting a new sequence: sample, sample_local, sample_global, sample_strelem')
    parser.add_argument("-k", "--kill", type=int, default=0, help='Timeout value of graph construction in seconds. (default: infinite)')
    parser.add_argument("-g", "--graphml", type=str, default=None, help='Write a graphml file with the given filename.')
    parser.add_argument("-c", "--csv", default=False, action='store_true', help='Write output as semi-colon csv file to stdout')
    parser.add_argument("-p", "--progress", default=False, action='store_true', help='Show progress of optimization')
    parser.add_argument("-d", "--debug", default=False, action='store_true', help='Show debug information of library')
    parser.add_argument("-r", "--range", type=str, default='6,20', help='Range of which a random number is generated to determine a spacer length')
    parser.add_argument("-t", "--three", type=int, default=10, help='Minimum length of the terminator stem, i.e. the number of nucleotides of the aptamer used to form a stem.')
    parser.add_argument("-u", "--ustretch", type=int, default=8, help='Length of the U stretch down stream of the terminator.')
    parser.add_argument("-x", "--xi", type=int, default=15, help='Energy that is added to the binding competent state to simulate ligand binding.')
    args = parser.parse_args()

    print("# Options: number={0:d}, jump={1:d}, exit={2:d}, strelem={3:d}, mode={4:}, nupack={5:}".format(args.number, args.jump, args.exit, args.strelem, args.mode, str(args.nupack)))
    rd.initialize_library(args.debug, args.kill)

    aptseq = ''
    aptstr = []
    spacersize = map(int, args.range.split(",")) #convert region to array of integers
    
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


    
    # main loop from zero to number of solutions
    for n in range(0, args.number):

        spacer = random.randint(spacersize[0],spacersize[1])
        spbp = int((spacer-4)/2)
        spss = (spacer-spbp*2)
        comp = len(aptseq) - random.randint(int(len(aptseq)/2),len(aptseq)-args.three+1) + 1
        
        # define structures
        structures = []
        constraint = ''
        start_sequence = ''
        print("# spacer length: " + str(spacer) + "\tlength of the complementary region: " + str(comp))
        ON = aptstr[0] + ("."*spacer) + ("."*comp) + ("."*args.ustretch)
        OFF=("."*(len(aptseq)-comp)) + ("("*comp) + ("("*spbp) + ("."*spss) + (")"*spbp) + (")"*comp)  + ("."*args.ustretch)
        if(spbp>3):
            OFF = ("."*(len(aptseq)-comp)) + ("("*comp) + ("."*(spbp-3)) + "(((" + ("."*spss) + ")))" + ("."*(spbp-3)) + (")"*comp)  + ("."*args.ustretch)
        
        SEQ = aptseq + ("N"*(spacer+comp))  + ("U"*args.ustretch)
        print("# Used input:\n# " + ON + "\n# " + OFF + "\n# " + SEQ)
        structures = [ON, OFF]
        constraint = SEQ
        
        # try to construct dependency graph, catch errors and timeouts
        dg = None
        construction_time = 0.0
        sample_time = 0.0
        
        # construct dependency graph with these structures
        try:
            start = time.clock()
            dg = rd.DependencyGraphMT(structures, constraint)
            construction_time = time.clock() - start
        except Exception as e:
            print( "Error: %s" % e , file=sys.stderr)
            
        # general DG values
        print("# " + "\n# ".join(structures) + "\n# " + constraint)
    
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
            # create a initial design object
            if (args.nupack):
                design = nupackDesign(structures, start_sequence)
            else:
                design = vrnaDesign(structures, start_sequence)
                
            # print header for csv file
            if (args.csv):
                print(";".join(["jump",
                                "exit",
                                "strelem",
                                "mode",
                                "score",
                                "num_mutations",
                                "construction_time",
                                "sample_time",
                                design.write_csv_header()] +
                               graph_properties.keys()))
    
            # reset the design object
            if (args.nupack):
                design = nupackDesign(structures, start_sequence)
            else:
                design = vrnaDesign(structures, start_sequence)
    
            #add variable
            design.aptstr=aptstr
            design.diff=args.xi
            design.aptsp_length=len(aptseq)+spacer
            
            start = time.clock()
            # do a complete sampling jump times
            (score, number_of_jumps) = classic_optimization(dg, design, objective_function=calculate_switch_objective, exit=args.jump, mode='sample', progress=args.progress)
            # now do the optimization based on the chosen mode
            try:
                (score, number_of_mutations) = classic_optimization(dg, design, objective_function=calculate_switch_objective, exit=args.exit, mode=args.mode, progress=args.progress)
            except ValueError as e:
                print (e.value)
                exit(1)
                # now do the optimization with mode strelem where we take structural elements and replace them a little
            number_of_strelem = 0
            if forgi_available:
                (score, number_of_strelem) = classic_optimization(dg, design, objective_function=calculate_switch_objective, exit=args.strelem, mode='sample_strelem', progress=args.progress)
            else:
                sys.stderr.write("-" * 60 + "\nWARNING: Strelem sampling not available!!!\nPlease install forgi https://github.com/pkerpedjiev/forgi\n" + "-" * 60 + "\n")
                sys.stderr.flush() 
                # sum up for a complete number of mutations
            number_of_mutations += number_of_jumps + number_of_strelem
            sample_time = time.clock() - start
    
            if (args.csv):
                print(args.jump,
                      args.exit,
                      args.strelem,
                      "\"" + args.mode + "\"",
                      score,
                      number_of_mutations,
                      construction_time,
                      sample_time,
                      design.write_csv(),
                      *graph_properties.values(), sep=";")
            else:
                print(design.write_out(score))
        else:
            print('# Construction time out reached!')

def calculate_switch_objective(design, weight=1):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    ((eos(1) - gibbs)+(aptsp_eos 0 gibbs_aptsp)/#structures) + 
                                    weight * abs((eos(1) - (eos(2)+xi)) +
    :param design: Design object containing the sequence and structures
    :param weight: To wheight the influence of the eos diffences
    '''

    aptsp_eos = RNA.energy_of_struct(design.sequence[:design.aptsp_length], design.structures[0][:design.aptsp_length])
    (aptsp_pf_str, aptsp_pf_energy) = RNA.pf_fold(design.sequence[:design.aptsp_length])
    
    objective_difference_part = abs(design.eos[0] - (design.eos[1]+design.diff))
    objective_ensemble = 0.5 * ((design.eos[1] - design.pf_energy) + (aptsp_eos - aptsp_pf_energy))
    (objective_bp_inboth, objective_bp_distance) = calculate_bp_distance(design.aptstr[0], design.mfe_structure)
    
    return (objective_ensemble + weight * (objective_difference_part) + objective_bp_inboth)


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
    
def create_bp_table(structure):
    bpo=[]
    bpt=[-1]*len(structure)
    for i, substr in enumerate(structure):
        if(substr=="("):
            bpo.append(i)
        elif(substr==")"):
            bpt[bpo.pop()] = i
    return bpt

if __name__ == "__main__":
    main()



