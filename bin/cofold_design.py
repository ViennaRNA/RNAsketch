from __future__ import print_function

try:
    from PyDesign import *
except ImportError, e:
    print(e.message)
    exit(1)

import RNAblueprint as rbp
import argparse
import sys
import time

def main():
    parser = argparse.ArgumentParser(description='Design a cofold device.')
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
    parser.add_argument("-r", "--reporter", type = str, default='CGTAAGGGCGAAGAGCTTTTTACCGGTGTTGTGCCTATTCTCGTAGAGTTAGATGGCGACGTTAAT', help='The coding sequence context, excluding the start codon that should be part of the sequence constraint. Default are the first 66 nucleotides of eGFP.')
    args = parser.parse_args()

    print("# Options: number={0:d}, exit={1:d}, mode={2:}, nupack={3:}".format(args.number, args.exit, args.mode, str(args.nupack)))
    rbp.initialize_library(args.debug, args.kill)
    # define structures
    structures = []
    fold_constraints = ['','']
    constraint = ''
    start_sequence = ''
    context = args.reporter
    if (args.input):
        data = ''
        for line in sys.stdin:
            data = data + '\n' + line
        (structures, constraint, start_sequence, fold_constraints) = read_input_additions(data)
    else:
        # RNAblueprint paper input of sRNA:5UTR design
        structures = [
            '........................................((((((((((((((((((((((((((&.............................)))))))))))))))))))))))))).........(((((((((((((((......))))))))))))))).....',
            '..................................................................&................................................................(((((((((((((((......))))))))))))))).....']
        fold_constraints = [
            '........................................((((((((((((((((((((((((((&.............................))))))))))))))))))))))))))..................................................',
            '.......................................xxxxxxxxxxxxxxxxxxxxxxxxxxx&............................xxxxxxxxxxxxxxxxxxxxxxxxxxxx<<<<<<<<.........................................']
        constraint = \
            'GGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAAGGAGNNNNNNNAUG&GGGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNUAGCAUAACCCCUUGGGGCCUCUAAACGGGUCUUGAGGGGUUUUUUG'
        context = \
            'CGTAAGGGCGAAGAGCTTTTTACCGGTGTTGTGCCTATTCTCGTAGAGTTAGATGGCGACGTTAAT'

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
        "[CU]{7,}",
        "GAAUUC",#EcoRI
        "UCUAGA",#XbaI
        "ACUAGU",#SpeI
        "CUGCAG" #PstI
    ]
    white_positions = [[m.start(), m.end()] for m in re.finditer('[AUGCT]+', constraint)]
    
    print( "# Parsed Input:\n# =============\n# structures: %s\n# fold constraints: %s\n# sequence constraints: %s\n# reporter context: %s\n# avoid motifs: %s\n# ignored positions: %s\n#" % (structures, fold_constraints, constraint, context,avoid_motifs,white_positions), file=sys.stdout)
    
    # try to construct dependency graph, catch errors and timeouts
    dg = None
    construction_time = 0.0
    sample_time = 0.0
        
    # construct dependency graph with these structures
    try:
        start = time.clock()
        dg = rbp.DependencyGraphMT(structures, constraint)
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
            print(";".join(["exit",
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
            design.foldconstraints = fold_constraints
            design.context = context
            #to evaluate binding site in standard output
            design.newState('binding', fold_constraints[0], constraint=fold_constraints[0])

            if (start_sequence):
                score=cofold_objective(design,printDetails=True)
                print(design.write_out(score))
            
            start = time.clock()
            # now do the optimization based on the chosen mode for args.exit iterations
            try:
                (score, number_of_mutations) = classic_optimization(dg, design, objective_function=cofold_objective, exit=args.exit, mode=args.mode, avoid_motifs=avoid_motifs, white_positions=white_positions, progress=args.progress)
            except ValueError as e:
                print (e.value)
                exit(1)
            # stop time counter    
            sample_time = time.clock() - start
            score=cofold_objective(design,printDetails=True)
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
    else:
        print('# Construction time out reached!')

def cofold_objective(design, weight1=1, weight2=1, weight3=1, printDetails=False):
    '''
    1 - [S AB ]/[A 0 ] + weight1 * P (RBS unpaired ) + weight2 * P(sRNA binding site unpaired) + weight3 * P(mRNA folds localy)
    1 - [AB] * e(-((Zsab/Zab')/(KT)))/ A0 + weight1 * P(RBS unpaired) + weight2 * P(sRNA binding site unpaired) + weight3 * P(mRNA folds localy)
    
    design object needs following inputexcept:
                raise 
    structure1 is the folded structure, in constraints we need only INTERmolecular pbs
    structure2 is the open state, in constraints we need the unpaired region (rbs) marked with xxxxxx and only mRNA!
    
    '''

    seqs = design.sequence.split('&')
    constr = design.foldconstraints[1].split('&')
    Ca0 = 1e-05

    # set RNA variables
    RNA.cvar.dangles = 2
    RNA.cvar.noLonelyPairs = 0
    RNA.cvar.cut_point = len(seqs[0])+1;
    
    # co_pf_fold = (structure, Gfe seq1, Gfe seq2, Gfe all INTER bp, Gfe all structs (dimers and monomers))
    # ab hetero dimers
    (x, ac, bc, fcab, cf) = RNA.co_pf_fold(seqs[0]+seqs[1]);
    # aa homo dimers
    (x, usel1, usel2, fcaa, usel3)= RNA.co_pf_fold(seqs[0]+seqs[0])
    # bb homo dimers
    RNA.cvar.cut_point = len(seqs[1])+1;
    (x, usel1, usel2, fcbb, usel3)= RNA.co_pf_fold(seqs[1]+seqs[1])

    # get concentration of Cab
    (Cab, Caa, Cbb, Ca, Cb)=RNA.get_concentrations(fcab, fcaa, fcbb, ac, bc, Ca0, 1e-03)
    RNA.cvar.cut_point = -1
    # save energy of all structures that build duplexes
    Eab = fcab
    # save gibbs free energy of mrna and sRNA
    Ea = ac
    Eb = bc

    # Get Energy of all structures that have our binding site Esab (constraint should look like: '....((((((.....&....))))))....')
    RNA.cvar.cut_point = len(seqs[0])+1;
    RNA.cvar.fold_constrained = 1
    (const_x, const_ac, const_bc, const_fcab, conqst_cf) = RNA.co_pf_fold(seqs[0]+seqs[1], remove_cuts(design.foldconstraints[0]))
    RNA.cvar.fold_constrained = 0
    RNA.cvar.cut_point = -1
    Esab = const_fcab
    
    # Get probability Psab of the duplex formed by the binding site
    Psab = Z_from_G(Esab-Eab)

    # Get concentration Csab of the duplex by multiplying Cab with the probability of Sab
    Csab = Cab * Psab

    # get probabilty that mRNA follows the given constraint (e.g RBS is unpaired '..xxxxxxxx..........')
    (con_str,con_energy)=getSaveConFold(seqs[0], constr[0],'pf_fold')
    PmRNAunpaired = Z_from_G(con_energy-Ea)
    
    # get probabilty that sRNA follows the given constraint (e.g sRNA binding side is unpaired '..xxxxxxxx..........')
    (con_str,con_energy)=getSaveConFold(seqs[1], constr[1],'pf_fold')
    PsRNAunpaired = Z_from_G(con_energy-Eb)

    # The designed 5'UTR sequence is extended at its 3' end by the
    # context, e.g. part of an reporter gene and the probability that
    # the mfe structure of the designed 5'UTR dominates the longer
    # construct
    (mRNA_structure,mRNA_energy) = getSaveConFold(seqs[0],constr[0],'fold')
    extendedSeq = (seqs[0] + design.context)
    extendedCon = (mRNA_structure.replace(".","x") + "."*len(design.context))
    (ext_pf_structure,ext_pf_energy) = getSaveConFold(extendedSeq,'','pf_fold')
    
    (con_str,con_energy)=getSaveConFold(extendedSeq, extendedCon,'pf_fold')
    PdesignDominates = Z_from_G(con_energy-ext_pf_energy)

    # calculate the probability that the 5'UTR structure and the
    # context structure folded individually dominate the structure
    # ensemble of the concatenated sequence
    (con_mRNA_pf_structure, con_mRNA_pf_energy) = getSaveConFold(seqs[0],constr[0],'pf_fold')
    (context_pf_structure, context_mRNA_pf_energy) = getSaveConFold(design.context,'','pf_fold')
    extendedSeq = (seqs[0] + design.context)
    (ext_pf_structure,ext_pf_energy) = getSaveConFold(extendedSeq,'','pf_fold')
    Pindividual = Z_from_G(con_mRNA_pf_energy+context_mRNA_pf_energy - ext_pf_energy)
    
    if(printDetails):
        print("Objective Details:\n=============\ncomplex concentration: %.2f P(5UTR_unpaired): %.2f P(sRNA_unpaired): %.2f P(mRNA_context): %.2f\n" % ((Csab / Ca0), PmRNAunpaired, PsRNAunpaired, Pindividual))
        print("1.0 - " + str(Csab) + "/" + str(Ca0) + " + " + str(weight1) + " * (1-" + str(PmRNAunpaired) + ") + " + str(weight2) + " * (1-" + str(PsRNAunpaired) + ") + " + str(weight3) + " * (1-" + str(Pindividual) + ")")
        print(str(1.0 - (Csab/Ca0)) + " + " + str(1-PmRNAunpaired) + " + " + str(1-PsRNAunpaired) + " + " + str(1-Pindividual))
        print(str(1.0 - Csab / Ca0 + weight1 * (1-PmRNAunpaired) + weight2 * (1-PsRNAunpaired) + weight3 * (1-Pindividual)))
    
    return 1.0 - Csab / Ca0 + weight1 * (1-PmRNAunpaired) + weight2 * (1-PsRNAunpaired) + weight3 * (1-Pindividual)

def getSaveConFold(sequence, constraint='', mode='fold'):
    # set RNA variables
    RNA.cvar.dangles = 2
    RNA.cvar.noLonelyPairs = 0
    RNA.cvar.fold_constrained = 1
    # pf_fold overwrites the constraint -> Bernard solved it!!! relys
    # on undefined behavior!!! we create a list from a string and then
    # join it into an presumably other string
    con = list(constraint)
    con = "".join(con)
    seq = list(sequence)
    seq = "".join(seq)
    
    #get free energy and structure
    (structure, energy) = ('',0.0)
    if(mode == 'fold'):
        (structure, energy) = RNA.fold(seq, con)
    elif(mode == 'pf_fold'):
        (structure, energy) = RNA.pf_fold(seq, con)
    else:
        print("Could not run getSaveConFold with mode" + mode, file=std.err)
        exit
    RNA.cvar.fold_constrained = 0
    # return free energy and structure
    return (structure, energy)


def getProbOfConstraintStructure(sequence, constraint, energy):
    # set RNA variables
    RNA.cvar.dangles = 2
    RNA.cvar.noLonelyPairs = 0
    RNA.cvar.fold_constrained = 1
    # pf_fold overwrites the constraint -> Bernard solved it!!! relys
    # on undefined behavior!!! we create a list from a string and then
    # join it into an presumably other string
    con = list(constraint)
    con = "".join(con)
    seq = list(sequence)
    seq = "".join(seq)
    #get gibbs free energy and ensemble structure
    (rna_x, rna_c) = RNA.pf_fold(seq, con)
    RNA.cvar.fold_constrained = 0
    # return the probability of the constraint being unpaired
    return Z_from_G(rna_c-energy)

def G_from_Z(Z):
    return - ((37.0 + 273.15)*1.98717)/1000.0 * math.log(Z)

def Z_from_G(G):
    return math.exp(- (G/ (((37.0 + 273.15)*1.98717)/1000.0)))

if __name__ == "__main__":
    main()


