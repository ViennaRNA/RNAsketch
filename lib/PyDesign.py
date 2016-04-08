#!/usr/bin/env python
'''
    PyDesign.py: A small module containing all the important helpers and 
    functions for all major RNAdesign operations.
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

import sys
import math
import re
import random
import collections

import RNAdesign as rd

'''
Global variable:
'''
vrna_available = True
nupack_available = True
forgi_available = True

try:
    import forgi.graph.bulge_graph as fgb
except ImportError, e:
    forgi_available = False

try:
    import RNA
except ImportError, e:
    vrna_available = False
    sys.stderr.write("-" * 60 + "\nWARNING: " + e.message + "!!!\n" + "-" * 60 + "\n")
    sys.stderr.flush()

try:
    import nupack
except ImportError, e:
    nupack_available = False
    sys.stderr.write("-" * 60 + "\nWARNING: " + e.message + "!!!\n" + "-" * 60 + "\n")
    sys.stderr.flush()

if not (vrna_available or nupack_available):
    raise ImportError("Neither ViennaRNA package nor Nupack found!")
    
class Design(object):
    '''
    Design contains all neccessary information to generate a RNA device
    or to return a solution. It is used as a container between the different
    functions.
    '''
    def __init__(self, structures, sequence=''):
        '''
        Construct a new Design object.

        :param sequence:
        :param structures:
        '''
        for struct in structures:
            create_bp_table(struct) #check for balanced brackets
        self.structures = structures
        self.sequence = sequence
        self._reset_all()
        self._temperatures = [37.0] * self.number_of_structures
        self._ligands = [None] * self.number_of_structures
        self._constraints = [None] * self.number_of_structures
    
    def _reset_sequence_dependent(self):
        self._eos = None
        self._pos = None
        self._eos_diff_mfe = None
        self._eos_reached_mfe = None
        self._mfe_structure = None
        self._mfe_energy = None
        self._pf_structure = None
        self._pf_energy = None
    
    def _reset_all(self):
        self._reset_sequence_dependent()
        self._number_of_structures = None
        self._length = None
        self._cut_points = None
        self._multifold = None
    
    @property
    def temperatures(self):
        return self._temperatures
    @temperatures.setter
    def temperatures(self, t):
        if isinstance(t, int):
            self._reset_sequence_dependent()
            self._temperatures = [t] * self.number_of_structures
        elif isinstance(t, list):
            if len(t) == self.number_of_structures:
                self._reset_sequence_dependent()
                self._temperatures = t
            else:
                raise TypeError('Temperature must be a list of doubles specifying the temperature for each structure')
        else:
            raise TypeError('Temperature must either be a list of doubles containing the temperature for every structure, or one integer.')
    
    @property
    def ligands(self):
        return self._ligands
    @ligands.setter
    def ligands(self, ligs):
        if isinstance(ligs, list) and len(ligs) == self.number_of_structures:
            for lig in ligs:
                if isinstance(lig, list) or not lig:
                    self._reset_sequence_dependent()
                    self._ligands = ligs
                else:
                    TypeError('A ligand list must either be None or must contain three values: sequence motif, struture motif and binding energy.')
        else:
            raise TypeError('Ligands must be a list of ligand lists for each structure.')
    
    @property
    def constraints(self):
        return self._constraints
    @constraints.setter
    def constraints(self, constraints):
        if isinstance(constraints, list) and len(constraints) == self.number_of_structures:
            for const in constraints:
                if isinstance(const, basestring) or not const:
                    if const:
                        create_bp_table(const) #check for balanced brackets
                    self._reset_sequence_dependent()
                    self._constraints = constraints
                else:
                    TypeError('A hard constraint must be a string containting only this characters: ().|x')
        else:
            raise TypeError('Constraints is a hard constraint string for each structure wraped in a list.')
    
    @property
    def classtype(self):
        return None
    
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, s):
        if isinstance(s, basestring) and re.match(re.compile("[AUGC\+\&]"), s):
            if len(s) != self.length:
                raise TypeError('Sequence must have the same length as the structural constraints')
            self._reset_sequence_dependent()
            self._sequence = s
        elif s == '':
            self._reset_sequence_dependent()
            self._sequence = None
        else:
            raise TypeError('Sequence must be a string containing a IUPAC RNA sequence')
    
    @property
    def structures(self):
        return self._structures
    @structures.setter
    def structures(self, s):
        if isinstance(s, list):
            length = len(s[0])
            for struct in s:
                if not (isinstance(struct, basestring) and re.match(re.compile("[\(\)\.\+\&]"), struct)):
                    raise TypeError('Structure be a string in dot-bracket notation')
                if length != len(struct):
                    raise TypeError('Structures must have equal length')
            self._reset_all()
            self._structures = s
        else:
            raise TypeError('Structures must be a list of dot-bracket strings')
    
    @property
    def eos(self):
        if not self._eos and self._sequence:
            self._eos = []
            for i, struct in enumerate(self.structures):
                self._eos.append(self._get_eos(self.sequence, struct, self.temperatures[i], self.ligands[i]))
        return self._eos
         
    @property
    def pos(self):
        if not self._pos and self._sequence:
            self._pos = []
            for i, eos in enumerate(self.eos):
                self._pos.append(math.exp((self.pf_energy[i]-eos) / self._get_KT(self.temperatures[i]) ))
        return self._pos
         
    @property
    def eos_diff_mfe(self):
        if not self._eos_diff_mfe and self._sequence:
            self._eos_diff_mfe = []
            for i, eos in enumerate(self.eos):
                self._eos_diff_mfe.append(eos - self.mfe_energy[i])
        return self._eos_diff_mfe
    
    @property
    def eos_reached_mfe(self):
        if not self._eos_reached_mfe and self._sequence:
            self._eos_reached_mfe = []
            for i, eos in enumerate(self.eos):
                if (eos == self.mfe_energy[i]):
                    self._eos_reached_mfe.append(1)
                else:
                    self._eos_reached_mfe.append(0)
        return self._eos_reached_mfe
    
    @property
    def mfe_energy(self):
        if not self._mfe_energy and self._sequence:
            self._calculate_mfe_energy_structure()
        return self._mfe_energy
    
    @property
    def mfe_structure(self):
        if not self._mfe_structure and self._sequence:
            self._calculate_mfe_energy_structure()
        return self._mfe_structure
    
    def _calculate_mfe_energy_structure(self):
        self._mfe_energy = []
        self._mfe_structure = []
        if self.temperatures[1:] == self.temperatures[:-1] and all(l is None for l in self.ligands) and all(c is None for c in self.constraints):
            (structure, energie) = self._get_fold(self.sequence, self.temperatures[0], self.ligands[0], self.constraints[0])
            self._mfe_energy = [energie] * self.number_of_structures
            self._mfe_structure = [structure] * self.number_of_structures
        else:
            for i, temperature in enumerate(self.temperatures):
                (structure, energie) = self._get_fold(self.sequence, temperature, self.ligands[i], self.constraints[i])
                self._mfe_energy.append(energie)
                self._mfe_structure.append(structure)
    
    @property
    def pf_energy(self):
        if not self._pf_energy and self._sequence:
            self._calculate_pf_energy_structure()
        return self._pf_energy
    
    @property
    def pf_structure(self):
        if not self._pf_structure and self._sequence:
            self._calculate_pf_energy_structure()
        return self._pf_structure
    
    def _calculate_pf_energy_structure(self):
        self._pf_energy = []
        self._pf_structure = []
        if self.temperatures[1:] == self.temperatures[:-1] and all(l is None for l in self.ligands) and all(c is None for c in self.constraints):
            (structure, energie) = self._get_pf_fold(self.sequence, self.temperatures[0], self.ligands[0], self.constraints[0])
            self._pf_energy = [energie] * self.number_of_structures
            self._pf_structure = [structure] * self.number_of_structures
        else:
            for i, temperature in enumerate(self.temperatures):
                (structure, energie) = self._get_pf_fold(self.sequence, temperature, self.ligands[i], self.constraints[i])
                self._pf_energy.append(energie)
                self._pf_structure.append(structure)
    
    def _get_KT(self, temperature):
        # KT = (betaScale*((temperature+K0)*GASCONST))/1000.0; /* in Kcal */
        return ((temperature + 273.15)*1.98717)/1000.0;
    
    def _get_eos(self, sequence, structure, temperature, ligand):
        raise NotImplementedError
    
    def _get_fold(self, sequence, temperature, ligand, constraint):
        raise NotImplementedError
    
    def _get_pf_fold(self, sequence, temperature, ligand, constraint):
        raise NotImplementedError
    
    @property
    def number_of_structures(self):
        if not self._number_of_structures:
            self._number_of_structures = len(self.structures)
        return self._number_of_structures
    
    @property
    def length(self):
        if not self._length:
            self._length = len(self.structures[0])
        return self._length
    
    @property
    def cut_points(self):
        if not self._cut_points:
            self._cut_points = []
            iterator = re.finditer(re.compile('\&|\+'), self.structures[0])
            for count, match in enumerate(iterator):
                self._cut_points.append(match.start()-count+1)
        return self._cut_points
    
    def _remove_cuts(self, input):
        return re.sub('[+&]', '', input)
    
    def _add_cuts(self, input):
        result = input
        for cut in self.cut_points:
            result = result[:cut-1] + '&' + result[cut-1:]
        return result
    
    @property
    def multifold(self):
        if not self._multifold:
            self._multifold = len(self.cut_points)
        return self._multifold
    
    def write_out(self, score=0):
        '''
        Generates a nice human readable version of all values of this design
        :param score: optimization score for this design
        :return: string containing a nicely formatted version of all design values
        '''
        result = '{0:}\t {1:5.2f}'.format(self.sequence, score)
        for i, struct in enumerate(self.structures):
            result += '\n{0:}\t{1:9.4f}\t{2:+9.4f}\t{3:9.4f}'.format(struct, self.eos[i], self.eos[i]-self.mfe_energy[i], self.pos[i])
            result += '\n{0:}\t{1:9.4f}'.format(self.mfe_structure[i], self.mfe_energy[i])
            result += '\n{0:}\t{1:9.4f}'.format(self.pf_structure[i], self.pf_energy[i])
        return result
    
    def write_csv(self, separator=';'):
        '''
        Generates a csv version of all values of this design separated by the given separator
        :param separator: separator for the values
        :return: string containing all values of this design separated by the given separator
        '''
        result = separator.join(map(str, ['\"' + self.sequence + '\"', self.length, self.number_of_structures ] + 
            self.mfe_energy +
            self.mfe_structure +
            self.pf_energy +
            self.pf_structure +
            self.eos + 
            self.eos_diff_mfe + 
            self.eos_reached_mfe + 
            self.pos))
        return result
    
    def write_csv_header(self, separator=';'):
        '''
        Generates a csv header for all values of this design separated by the given separator
        :param separator: separator for the values
        :return: string containing a csv header for this design separated by the given separator
        '''
        result = separator.join(['sequence', 'seq_length', 'number_of_structures'])
        strings = ['mfe_energy_', 'mfe_structure_', 'pf_energy_', 'pf_structure_', 'eos_', 'diff_eos_mfe_', 'mfe_reached_', 'prob_']
        for s in strings:
            for i in range(0, self.number_of_structures):
                result += separator + s + str(i)
        return result
    
    def _create_header(self, separator, strings):
        result = ''
        
        return result

if vrna_available:
    class vrnaDesign(Design):
        @property
        def classtype(self):
            return 'vrnaDesign'
        
        def _change_cuts(self, input):
            return re.sub('[+]', '&', input)
        
        def _get_eos(self, sequence, structure, temperature, ligand=None):
            if self.multifold > 1:
                raise NotImplementedError
            RNA.cvar.temperature = temperature
            fc = RNA.fold_compound(self._change_cuts(sequence), None, RNA.VRNA_OPTION_MFE | RNA.VRNA_OPTION_EVAL_ONLY)
            if ligand:
                fc.sc_add_hi_motif(ligand[0], ligand[1], ligand[2])
            return fc.eos(self._remove_cuts(structure))
    
        def _get_fold(self, sequence, temperature, ligand=None, constraint=None):
            RNA.cvar.temperature = temperature
            fc = RNA.fold_compound(self._change_cuts(sequence), None, RNA.VRNA_OPTION_MFE)
            if ligand:
                fc.sc_add_hi_motif(ligand[0], ligand[1], ligand[2])
            if constraint:
                fc.hc_add_db(self._remove_cuts(constraint))
            if self.multifold == 0:
                (structure, energie) = fc.mfe()
            if self.multifold == 1:
                (structure, energie) = fc.mfe_dimer()
                structure = self._add_cuts(structure)
            if self.multifold > 1:
                raise NotImplementedError
            return (structure, energie)
    
        def _get_pf_fold(self, sequence, temperature, ligand=None, constraint=None):
            RNA.cvar.temperature = temperature
            fc = RNA.fold_compound(self._change_cuts(sequence), None, RNA.VRNA_OPTION_PF)
            if ligand:
                fc.sc_add_hi_motif(ligand[0], ligand[1], ligand[2])
            if constraint:
                fc.hc_add_db(self._remove_cuts(constraint))
            if self.multifold == 0:
                (structure, energie) = fc.pf()
            if self.multifold == 1:
                (structure, energie) = fc.pf_dimer()
                structure = self._add_cuts(structure)
            elif self.multifold > 1:
                raise NotImplementedError
            return (structure, energie)

if nupack_available:
    class nupackDesign(Design):
        @property
        def classtype(self):
            return 'nupackDesign'
        
        def _change_cuts(self, input):
            return re.sub('[&]', '+', input)
    
        def _get_eos(self, sequence, structure, temperature, ligand=None):            
            #TODO nupack.energy can not handle unconnected cofold structures
            return nupack.energy([self._change_cuts(sequence)], self._change_cuts(structure), material = 'rna', pseudo = True, T = temperature)
    
        def _get_fold(self, sequence, temperature, ligand=None, constraint=None):
            nupack_mfe = nupack.mfe([self._change_cuts(sequence)], material = 'rna', pseudo = True, T = temperature) # if str, 0, no error

            pattern = re.compile('(\[\(\')|(\',)|(\'\)\])')
            temp_mfe = pattern.sub('', "%s" %nupack_mfe)
            temp_mfe = temp_mfe.replace("'", "")
            mfe_list = temp_mfe.split()
    
            mfe_struct = mfe_list[0]
            mfe_energy = float(mfe_list[1])
            return mfe_struct, mfe_energy
    
        def _get_pf_fold(self, sequence, temperature, ligand=None, constraint=None):
            # Nupack doesn't return ensemble structure
            return re.sub('[^\+]', '?', self._change_cuts(sequence)), nupack.pfunc([sequence], material = 'rna', pseudo = True, T = temperature)

def read_inp_file(filename):
    '''
    Reads a file in *.inp format and returns all neccessary information
    :param filename: Filename of the file to read
    :return structures: List of structures in dot-bracket notation
    :return constraint: Sequence constraint
    '''
    
    with open(filename) as f:
        data = f.read()
        (structures, constraint, sequence) = read_input(data)
    
    constraint = constraint.replace(" ", "N")
    return structures, constraint, sequence

def read_input(content):
    '''
    Reads some input and returns all neccessary information in the right container.
    Input is a string, lines separated by linebreaks. Content might be structures, 
    a sequence constraint and a start sequence
    :param filename: Filename of the file to read
    :return structures: List of structures in dot-bracket notation
    :return constraint: Sequence constraint
    :return sequence: Start sequence
    '''
    structures = []
    constraint = ''
    sequence = ''
    
    lines = content.split("\n")
    for line in lines:
        if re.match(re.compile("^[\(\)\.\{\}\[\]\<\>\+\&]+$"), line, flags=0):
            structures.append(line.rstrip('\n'))
        elif re.match(re.compile("^[\ ACGTUWSMKRYBDHVN\&\+]+$"), line, flags=0):
            line = line.replace(" ", "N")
            if re.match(re.compile("^[ACGTU\&\+]+$"), line, flags=0) and sequence == '':
                sequence = line.rstrip('\n')
            elif constraint == '':
                constraint = line.rstrip('\n')
            else:
                raise ValueError('Too many constraints or start sequences: ' + line)
    
    checklength = len(structures[0])
    if constraint != '' and len(constraint) != checklength:
        raise ValueError('Structures and the sequence constraint must have the same length!')
    elif sequence != '' and len(sequence) != checklength:
        raise ValueError('Structures and the start sequence must have the same length!')
    
    for s in structures:
        if len(s) != checklength:
            raise ValueError('Structures must all have the same length!')
    
    return structures, constraint, sequence

def create_bp_table(structure):
    '''
    Takes a structure in dot bracket notation and returns a base pair table.
    Unpaired positions are -1, otherwise the index of the adjacent bracket is listed
    :param structure: string with dot-bracket notation of the strcture
    :return bpt: base pair table
    '''
    bpo=[]
    bpt=[-1]*len(structure)
    for i, substr in enumerate(structure):
        if(substr=="("):
            bpo.append(i)
        elif(substr==")"):
            try:
                bpt[bpo.pop()] = i
            except:
                raise ValueError('Unbalanced brackets: too few opening brackets')
    if len(bpo) > 0:
        raise LogicError('Unbalanced brackets: too few closing brackets')
    return bpt

def get_graph_properties(dg):
    '''
    Takes a RNAdesign DependencyGraph Object and constructs a dicionary with all the
    calculated properties.
    :param dg: RNAdesign DependencyGraph object
    :return properties: Dictionary containing all the graph properties
    '''
    properties = {}
    special_ratios = []
    max_specials = 0
    max_component_vertices = 0
    
    properties['num_cc'] = dg.number_of_connected_components()
    properties['nos'] = dg.number_of_sequences()
    
    for cc in range(0, properties['num_cc']):
        cv = len(dg.component_vertices(cc))
        sv = len(dg.special_vertices(cc))
        special_ratios.append(float(sv)/float(cv))
        if (max_specials < sv):
            max_specials = sv
        if (max_component_vertices < cv):
            max_component_vertices = cv
    
    properties['max_specials'] = max_specials
    properties['max_component_vertices'] = max_component_vertices
    properties['max_special_ratio'] = max(special_ratios)
    properties['mean_special_ratio'] = sum(special_ratios) / len(special_ratios)
        
    return properties

def calculate_objective(design, weight=0.5):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    (eos(1)+eos(2)+eos(3) - 3 * gibbs) / number_of_structures +
        weight * (eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2) * 2 / (number_of_structures * (number_of_structures-1))
    :param design: Design object containing the sequence and structures
    :param weight: To wheight the influence of the eos diffences
    '''
    return calculate_objective_1(design) + weight * calculate_objective_2(design)

def calculate_objective_1(design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    (eos(1)+eos(2)+eos(3) - 3 * gibbs) / number_of_structures
    :param design: Design object containing the sequence and structures
    '''
    return (sum(design.eos) - sum(design.pf_energy)) / design.number_of_structures

def calculate_objective_2(design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    (eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2) * 2 / (number_of_structures * (number_of_structures-1))
    :param design: Design object containing the sequence and structures
    '''
    objective_difference_part = 0
    for i, eos1 in enumerate(design.eos):
        for eos2 in design.eos[i+1:]:
            objective_difference_part += math.fabs(eos1 - eos2)
    
    return objective_difference_part * 2 / (design.number_of_structures * (design.number_of_structures-1))

def sample_sequence(dg, design, mode, sample_steps=1, avoid=[]):
    '''
    This function samples a sequence with the given mode from the dependency graph object
    and writes it into the design object
    :param dg: RNAdesign dependency graph object
    :param design: design object
    :param mode: mode how to sample, this is a string
    :param sample_steps: count how many times to do the sample operation
    :param return: mut_nos is the solution space we drew from
    :param return: sample_count is how many times we sampled a solution from the dependency graph object (important for revert later)
    '''
    # count how many samples we did to be able to revert this later
    sample_count = 0
    # remember the solution space we drew from
    mut_nos = 1
    dg.set_history_size(sample_count + 100)
    
    # sample a new sequence
    # if random choice is requested pick something new
    if mode == "random":
        modes = ['sample','sample_global','sample_local', 'sample_strelem']
        mode = random.choice(modes)
    if sample_steps == 0:
        sample_steps = random.randrange(1, dg.number_of_connected_components())

    if mode == 'sample':
        mut_nos = dg.sample()
        sample_count += 1
    elif mode == 'sample_global':
        for c in _sample_connected_components(dg, sample_steps):
            mut_nos *= dg.sample_global(c)
            sample_count += 1
    elif mode == 'sample_local':
        # TODO this local sampling is unfair this way and mut_nos is not calculated correctly!
        for _ in range(0, sample_steps):
            mut_nos *= dg.sample_local()
            sample_count += 1
    elif mode == 'sample_strelem':
        if forgi_available:
            # sample new sequences for structural elements
            struct = design._remove_cuts(random.choice(design.structures))
            bg = fgb.BulgeGraph(dotbracket_str=struct)
            for s in bg.random_subgraph(sample_steps):
                try:
                    mut_nos *= dg.sample(bg.defines[s][0]-1, bg.defines[s][1]-1)
                    sample_count += 1
                    if s[0] == 'i':
                        mut_nos *= dg.sample(bg.defines[s][2]-1, bg.defines[s][3]-1)
                        sample_count += 1
                except IndexError:
                    pass
        else:
            raise ImportError("Forgi Library not available!")
    else:
        raise ValueError("Wrong mode argument: " + mode + "\n")
    
    # assign sequence to design and return values
    design.sequence = dg.get_sequence()
    return (mut_nos, sample_count)

def classic_optimization(dg, design, objective_function=calculate_objective, exit=1000, mode='sample', progress=False):
    '''
    Takes a Design object and does a classic optimization of this sequence.
    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_functions: array of functions which takes a design object and returns a score for evaluation
    :param exit: Number of unsuccessful new sequences before exiting the optimization
    :param mode: String defining the sampling mode: sample, sample_global, sample_local
    :param progress: Whether or not to print the progress to the console
    :param return: Optimization score reached for the final sequence
    "param return: Number of samples neccessary to reach this result
    '''
    # if the design has no sequence yet, sample one from scratch
    if not design.sequence:
        dg.sample()
        design.sequence = dg.get_sequence()
    else:
        dg.set_sequence(design.sequence)
    
    score = objective_function(design)
    # count for exit condition
    count = 0
    # remember how may mutations were done
    number_of_samples = 0
    
    # main optimization loop 
    while True:
        # count up the mutations
        number_of_samples += 1
        # sample a new sequence
        (mut_nos, sample_count) = _sample_sequence(dg, design, mode)
        
        # write progress
        if progress:
            sys.stderr.write("\rMutate: {0:7.0f}/{1:5.0f} | Score: {2:5.2f} | NOS: {3:.5e}".format(number_of_samples, count, score, mut_nos) + " " * 20)
            sys.stderr.flush()
        
        this_score = objective_function(design)
        # evaluate
        if (this_score < score):
            score = this_score
            count = 0
        else:
            dg.revert_sequence(sample_count)
            design.sequence = dg.get_sequence()
            count += 1
            if count > exit:
                break
    
    # clear the console
    if (progress):
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()
    
    # finally return the result
    return score, number_of_samples

def constraint_generation_optimization(dg, design, objective_function=calculate_objective, exit=1000, mode='sample', num_neg_constraints=100, max_eos_diff=0, progress=False):
    '''
    Takes a Design object and does a constraint generation optimization of this sequence.
    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_functions: array of functions which takes a design object and returns a score for evaluation
    :param exit: Number of unsuccessful new sequences before exiting the optimization
    :param mode: String defining the sampling mode: sample, sample_global, sample_local
    :param num_neg_constraints: Maximal number of negative constraints to accumulate during the optimization process
    :param max_eos_diff: Maximal difference between eos of the negative and positive constraints
    :param progress: Whether or not to print the progress to the console
    :param return: Optimization score reached for the final sequence
    "param return: Number of samples neccessary to reach this result
    '''
    dg.set_history_size(100)
    neg_constraints = collections.deque(maxlen=num_neg_constraints)
    
    # if the design has no sequence yet, sample one from scratch
    if not design.sequence:
        dg.sample()
        design.sequence = dg.get_sequence()
    else:
        dg.set_sequence(design.sequence)
    
    score = objective_function(design)
    # count for exit condition
    count = 0
    # remember how may mutations were done
    number_of_samples = 0
    
    # main optimization loop
    while True:
        # constraint generation loop
        while True:
            # count up the mutations
            number_of_samples += 1
            # sample a new sequence
            (mut_nos, sample_count) = _sample_sequence(dg, design, mode)
            
            # write progress
            if progress:
                sys.stderr.write("\rMutate: {0:7.0f}/{1:5.0f} | EOS-Diff: {2:4.2f} | Score: {3:5.2f} | NOS: {4:.5e}".format(number_of_samples, count, max_eos_diff, score, mut_nos))
                sys.stderr.flush()
            # boolean if it is perfect already
            perfect = True
            # evaluate the constraints
            for negc in reversed(neg_constraints):
                # test if the newly sampled sequence is compatible to the neg constraint, if not -> Perfect!
                if rd.sequence_structure_compatible(design.sequence, [negc]):
                    if design.classtype == 'vrnaDesign':
                        neg_eos = RNA.energy_of_struct(design.sequence, negc)
                    elif design.classtype == 'nupackDesign':
                        neg_eos = nupack.energy([design.sequence], negc, material = 'rna', pseudo = True)
                    else:
                        raise ValueError('Could not figure out the classtype of the Design object.')
                    # test if the newly sampled sequence eos for pos constraints is lower than
                    # the eos for all negative constraints, if not -> Perfect!
                    for k in range(0, design.number_of_structures):
                        if float(neg_eos) - float(design.eos[k]) < max_eos_diff:
                            # this is no better solution, revert!
                            perfect = False
                            dg.revert_sequence(sample_count)
                            design.sequence = dg.get_sequence()
                            break
                # if this is no perfect solution, stop evaluating and sample a new one
                if not perfect:
                    break
            # if solution is perfect, stop the optimization and go down to score calculation
            if perfect:
                break
        
        # count this as a solution to analyse
        count += 1
        # calculate objective
        this_score = objective_function(design)

        if (this_score < score):
            score = this_score
            # reset values
            count = 0
        else:
            dg.revert_sequence(sample_count)
            design.sequence = dg.get_sequence()
        # else if current mfe is not in negative constraints, add to it
        for mfe_str in design.mfe_structure:
            if mfe_str not in design.structures:
                if mfe_str not in neg_constraints:
                    neg_constraints.append(mfe_str)
                    #print('\n'+'\n'.join(neg_constraints))

        # exit condition
        if count > exit:
            break
        
    # clear the console
    if (progress):
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()
    # finally return the result
    return score, number_of_samples

def _sample_connected_components(dg, amount=1):
    '''
    This function samples several connected component weighted by their number of solutions.
    We need this function to draw from the set of CCs without getting the same CC twice.
    Therefore the complete number of solutions if we sample these CCs new, is the product of the
    number of solutions for each CC.
    :param dg: Dependency Graph object from the RNAdesig library
    :param amount: number of connected components to sample
    :return: list of connected component IDs which can be used for example for: dg.sample_global(ID)
    '''
    result = []
    noslist = {}
    
    if amount > dg.number_of_connected_components():
        amount = dg.number_of_connected_components()
    
    for c in range(0, dg.number_of_connected_components()):
        noslist[c] = dg.number_of_sequences(c)
    
    for _ in range(0, amount):
        rand = random.randint(0, sum(noslist.values())-1)
        keys = []
        for c in noslist.keys():
            keys.append(c)
            if rand < sum([noslist[key] for key in keys]):
                result.append(c)
                del noslist[c]
                break
    return result

def sample_count_unique_solutions(solution_space_size, sample_size):
    '''
    Calculates the expectancy value of how many time it is necessary to draw a 
    solution to gain a unique set of solutions with the given sample size
    :param solution_space_size: The size of the complete solutions space to draw from
    :param sample_size: The size of the requested unique set
    :return: Expectancy value of how many times to draw from the solution space to gain this unique set
    '''
    sum = float(0)
    for k in xrange(int(solution_space_size) - int(sample_size) + 1, int(solution_space_size) + 1):
        sum += float(solution_space_size) / float(k)
    return float(sum)

