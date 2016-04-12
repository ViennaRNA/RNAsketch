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
    parser = argparse.ArgumentParser(description='Design a cofold device.')
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
    args = parser.parse_args()

    print("# Options: number={0:d}, jump={1:d}, exit={2:d}, strelem={3:d}, mode={4:}, nupack={5:}".format(args.number, args.jump, args.exit, args.strelem, args.mode, str(args.nupack)))
    rd.initialize_library(args.debug, args.kill)
    # define structures
    structures = []
    fold_constraints = []
    constraint = ''
    start_sequence = ''
    
    if (args.input):
        data = ''
        for line in sys.stdin:
            data = data + '\n' + line
        (structures, constraint, start_sequence) = read_input(data)
    elif (args.file is not None):
        print("# Input File: {0:}".format(args.file))
        (structures, constraint, start_sequence) = read_inp_file(args.file)
    else:
        structures = [
            '........................................((((((((((((((((((((((((((((((.......................................&...((((((((((.....)))))))))).)))))))))))))))))))))))))))))).........(((((((((((((((......))))))))))))))).....',
            '......((((((((((((((......)))))))))))))).....................................................................&...((((((((((.....))))))))))........................................(((((((((((((((......))))))))))))))).....']
        fold_constraints = [
            '........................................((((((((((((((((((((((((((((((.......................................&.............................))))))))))))))))))))))))))))))..................................................',
            '..............................................xxxxxxxxxxxxxxxxxxxxxxxx.......................................&.............................................................................................................']
        constraint = \
            'GGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGGAGNNNNNNNAUGAACAGCAGCAACCUGGCGGCAGCGCAAAAGAUGCGUAAA&GGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAGCAUAACCCCUUGGGGCCUCUAAACGGGUCUUGAGGGGUUUUUUG'
    
    #sequence motifs to avoid
    avoid_motifs = [
        "[A]{4,}",
        "[C]{4,}",
        "[G]{4,}",
        "[U]{4,}",
        "[GU]{7,}",
        "[AC]{7,}",
        "[AG]{7,}",
        "[GC]{7,}",
        "[AU]{7,}",
        "[CU]{7,}"
    ]
    
    white_positions = [[m.start(), m.end()] for m in re.finditer('[AUGCT]+', constraint)]
    print(white_positions)
    
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

        # main loop from zero to number of solutions
        for n in range(0, args.number):
            # reset the design object
            if (args.nupack):
                design = nupackDesign(structures, start_sequence)
            else:
                design = vrnaDesign(structures, start_sequence)
            
            # set fold constraints
            design.constraints = fold_constraints
            
            start = time.clock()
            # do a complete sampling jump times
            (score, number_of_jumps) = classic_optimization(dg, design, objective_function=cofold_objective, exit=args.jump, mode='sample', avoid_motifs=avoid_motifs, white_positions=white_positions, progress=args.progress)
            # now do the optimization based on the chose mode
            try:
                (score, number_of_mutations) = classic_optimization(dg, design, objective_function=cofold_objective, exit=args.exit, mode=args.mode, avoid_motifs=avoid_motifs, white_positions=white_positions, progress=args.progress)
            except ValueError as e:
                print (e.value)
                exit(1)
            # now do the optimization with mode strelem where we take structural elements and replace them a little
            number_of_strelem = 0
            if forgi_available:
                (score, number_of_strelem) = classic_optimization(dg, design, objective_function=cofold_objective, exit=args.strelem, mode='sample_strelem', avoid_motifs=avoid_motifs, white_positions=white_positions, progress=args.progress)
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
                print(design.sequence)
                print(design.write_out(score))
    else:
        print('# Construction time out reached!')

def cofold_objective(design, weight=1):
    '''
    1 - [S AB ]/[A 0 ] + weight * P (RBS unpaired )
    1 - [AB] * e(-((Zsab/Zab')/(KT)))/ A0 + weight * P(RBS unpaired)
    
    design object needs following inputexcept:
                raise 
    structure1 is the folded structure, in constraints we need only INTERmolecular pbs
    structure2 is the open state, in constraints we need the unpaired region (rbs) marked with xxxxxx and only mRNA!
    
    '''
    #print('seq = ' + design.sequence)
    seqs = design.sequence.split('&')
    constr = design.constraints[1].split('&')
    Ca0 = 1e-05
    
    # get concentration of ab: Cab
    # ab
    RNA.cvar.dangles = 2
    RNA.cvar.noLonelyPairs = 0

    RNA.cvar.cut_point = len(seqs[0])+1;
    # structure, Gfe seq1, Gfe seq2, Gfe all INTER bp, Gfe all structs (dimers and monomers)
    (x, ac, bc, fcab, cf) = RNA.co_pf_fold(seqs[0]+seqs[1]);
    #print(x)
    # aa
    (x, usel1, usel2, fcaa, usel3)= RNA.co_pf_fold(seqs[0]+seqs[0])
    # bb
    RNA.cvar.cut_point = len(seqs[1])+1;
    (x, usel1, usel2, fcbb, usel3)= RNA.co_pf_fold(seqs[1]+seqs[1])
    
    (Cab, Caa, Cbb, Ca, Cb)=RNA.get_concentrations(fcab, fcaa, fcbb, ac, bc, Ca0, 1e-03);
    RNA.cvar.cut_point = -1
    #print('Cab = ' + str(Cab))
    # save energy of all structures that build duplexes
    Eab = fcab
    #print('Eab = ' + str(Eab))
    # save gibbs free energy of mrna
    Ea = ac
    #print('Ea = ' + str(Ea))
    
    # get Esab  (constraint should look like: '....((((((.....&....))))))....')
    RNA.cvar.cut_point = len(seqs[0])+1;
    RNA.cvar.fold_constrained = 1
    (const_x, const_ac, const_bc, const_fcab, conqst_cf) = RNA.co_pf_fold(seqs[0]+seqs[1], design._remove_cuts(design.constraints[0]));
    RNA.cvar.fold_constrained = 0
    #print(const_x)
    RNA.cvar.cut_point = -1
    # Energy of all structures that have our binding site
    Esab = const_fcab
    #print('Esab = ' + str(Esab))
    
    # get Psab
    Psab = Z_from_G(Esab-Eab)
    #print('Psab = ' + str(Psab))
    # get concentration of Sab by multiplying Cab with the probability of Sab
    Csab = Cab * Psab
    #print('Csab = ' + str(Csab))
    # get Prbs unpaired (constraint should look like: '..xxxxxxxx..........')
    RNA.cvar.fold_constrained = 1
    (mrna_x, mrna_c) = RNA.pf_fold(seqs[0], constr[0]); #TODO constraints only length of seqs[0]
    RNA.cvar.fold_constrained = 0
    # save the energy of all structures with the RBS being unpaired
    Erbsunpaired = mrna_c
    #print('Erbsunpaired = ' + str(Erbsunpaired))
    
    # get the probability of the RBS being unpaired
    Prbsunpaired =  Z_from_G(Erbsunpaired-Ea)
    #print('Prbsunpaired = ' + str(Prbs))
    # calculate the objective and return
    #print('score = ' + str(1.0 - Csab / Ca0 + weight * Prbsunpaired))
    return 1.0 - Csab / Ca0 + weight * (1-Prbsunpaired)

def G_from_Z(Z):
    return - ((37.0 + 273.15)*1.98717)/1000.0 * math.log(Z)

def Z_from_G(G):
    return math.exp(- (G/ (((37.0 + 273.15)*1.98717)/1000.0)))

if __name__ == "__main__":
    main()


