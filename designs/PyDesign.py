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
# KT = (betaScale*((temperature+K0)*GASCONST))/1000.0; /* in Kcal */
KT = ((37+273.15)*1.98717)/1000.0;
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
        self.structures = structures
        self.sequence = sequence
        self._reset_values()
    
    def _reset_values(self):
        self._eos = None
        self._pos = None
        self._eos_diff_mfe = None
        self._eos_reached_mfe = None
        self._mfe_structure = None
        self._mfe_energy = None
        self._pf_structure = None
        self._pf_energy = None
        self._number_of_structures = None
        self._length = None
    
    @property
    def classtype(self):
        return None
    
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, s):
        if isinstance(s, basestring) and re.match(re.compile("[AUGC]"), s):
            if len(s) != self.length:
                raise TypeError('Sequence must have the same length as the structural constraints')
            self._reset_values()
            self._sequence = s
        elif s == '':
            self._reset_values()
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
                if not (isinstance(struct, basestring) and re.match(re.compile("[\(\)\.]"), struct)):
                    raise TypeError('Structure be a string in dot-bracket notation')
                if length != len(struct):
                    raise TypeError('Structures must have equal length')
            self._reset_values()
            self._structures = s
        else:
            raise TypeError('Structures must be a list of dot-bracket strings')
    
    @property
    def eos(self):
        if not self._eos and self._sequence:
            self._eos = []
            for struct in self.structures:
                self._eos.append(self._get_eos(self.sequence, struct))
        return self._eos
         
    @property
    def pos(self):
        if not self._pos and self._sequence:
            self._pos = []
            for eos in self.eos:
                self._pos.append(math.exp((self.pf_energy-eos) / KT ))
        return self._pos
         
    @property
    def eos_diff_mfe(self):
        if not self._eos_diff_mfe and self._sequence:
            self._eos_diff_mfe = []
            for eos in self.eos:
                self._eos_diff_mfe.append(eos - self.mfe_energy)
        return self._eos_diff_mfe
    
    @property
    def eos_reached_mfe(self):
        if not self._eos_reached_mfe and self._sequence:
            self._eos_reached_mfe = []
            for eos in self.eos:
                if (eos == self.mfe_energy):
                    self._eos_reached_mfe.append(1)
                else:
                    self._eos_reached_mfe.append(0)
        return self._eos_reached_mfe
    
    @property
    def mfe_energy(self):
        if not self._mfe_energy and self._sequence:
            (self._mfe_structure, self._mfe_energy) = self._get_fold(self.sequence)
        return self._mfe_energy
    
    @property
    def mfe_structure(self):
        if not self._mfe_structure and self._sequence:
            (self._mfe_structure, self._mfe_energy) = self._get_fold(self.sequence)
        return self._mfe_structure
    
    @property
    def pf_energy(self):
        if not self._pf_energy and self._sequence:
            (self._pf_structure, self._pf_energy) = self._get_pf_fold(self.sequence)
        return self._pf_energy
    
    @property
    def pf_structure(self):
        if not self._pf_structure and self._sequence:
            (self._pf_structure, self._pf_energy) = self._get_pf_fold(self.sequence)
        return self._pf_structure
    
    def _get_eos(self, sequence, structure):
        raise NotImplementedError
    
    def _get_fold(self, sequence):
        raise NotImplementedError
    
    def _get_pf_fold(self):
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
    
    def write_out(self, score=0):
        '''
        Generates a nice human readable version of all values of this design
        :param score: optimization score for this design
        :return: string containing a nicely formatted version of all design values
        '''
        result = '{0:}\t{1:9.4f}'.format(self.sequence, score)
        for i, struct in enumerate(self.structures):
            result += '\n{0:}\t{1:9.4f}\t{2:+9.4f}\t{3:9.4f}'.format(struct, self.eos[i], self.eos[i]-self.mfe_energy, self.pos[i])
        result += '\n{0:}\t{1:9.4f}'.format(self.mfe_structure, self.mfe_energy)
        result += '\n{0:}\t{1:9.4f}'.format(self.pf_structure, self.pf_energy)
        return result
    
    def write_csv(self, separator=';'):
        '''
        Generates a csv version of all values of this design separated by the given separator
        :param separator: separator for the values
        :return: string containing all values of this design separated by the given separator
        '''
        result = separator.join(map(str, ['\"' + self.sequence + '\"', self.length, self.number_of_structures, self.mfe_energy, 
            '\"'+ self.mfe_structure + '\"', self.pf_energy, '\"' + self.pf_structure + '\"'] + 
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
        result = separator.join(['sequence', 'seq_length', 'number_of_structures', 'mfe_energy', 'mfe_structure', 'pf_energy', 'pf_structure'])
        eos_string = ''
        pos_string = ''
        eos_diff_mfe_string = ''
        eos_reached_mfe_string = ''
        for s in range(0, self.number_of_structures):
            eos_string += separator + "eos_" + str(s)
            eos_diff_mfe_string += separator + "diff_eos_mfe_" + str(s)
            eos_reached_mfe_string += separator + "mfe_reached_" + str(s)
            pos_string += separator + "prob_" + str(s)
        result += eos_string + eos_diff_mfe_string + eos_reached_mfe_string + pos_string
        return result

if vrna_available:
    class vrnaDesign(Design):
        @property
        def classtype(self):
            return 'vrnaDesign'
    
        def _get_eos(self, sequence, structure):
            return RNA.energy_of_struct(sequence, structure)
    
        def _get_fold(self, sequence):
            return RNA.fold(self.sequence)
    
        def _get_pf_fold(self, sequence):
            return RNA.pf_fold(self.sequence)

if vrna_available:
    class nupackDesign(Design):
        @property
        def classtype(self):
            return 'nupackDesign'
    
        def _get_eos(self, sequence, structure):
            return nupack.energy([sequence], structure, material = 'rna', pseudo = True)
    
        def _get_fold(self, sequence):
            nupack_mfe = nupack.mfe([sequence], material = 'rna', pseudo = True) # if str, 0, no error
            pattern = re.compile('(\[\(\')|(\',)|(\'\)\])')
            temp_mfe = pattern.sub('', "%s" %nupack_mfe)
            temp_mfe = temp_mfe.replace("'", "")
            mfe_list = temp_mfe.split()
    
            mfe_struct = mfe_list[0]
            mfe_energy = float(mfe_list[1])
            return mfe_struct, mfe_energy
    
        def _get_pf_fold(self, sequence):
            # Nupack doesn't return ensemble structure
            return '?' * len(sequence), nupack.pfunc([sequence], material = 'rna', pseudo = True)

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
        if re.match(re.compile("^[\(\)\.\{\}\[\]\<\>]+$"), line, flags=0):
            structures.append(line.rstrip('\n'))
        elif re.match(re.compile("^[\ ACGTUWSMKRYBDHVN]+$"), line, flags=0):
            line = line.replace(" ", "N")
            if re.match(re.compile("^[ACGTU]+$"), line, flags=0) and sequence == '':
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

def calculate_objective(design, weight=1):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    eos(1)+eos(2)+eos(3) - 3 * gibbs + 
                                    weight * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2) / (3!/(3-2)!*2)
    :param design: Design object containing the sequence and structures
    :param weight: To wheight the influence of the eos diffences
    '''
    objective_difference_part = 0
    for i, eos1 in enumerate(design.eos):
        for eos2 in design.eos[i+1:]:
            objective_difference_part += math.fabs(eos1 - eos2)
    
    combination_count = 1
    if (design.number_of_structures != 1):
        combination_count = math.factorial(design.number_of_structures) / (math.factorial(design.number_of_structures-2)*2)
    
    return sum(design.eos) - design.number_of_structures * design.pf_energy + weight * (objective_difference_part / combination_count)

def _sample_sequence(dg, design, mode, sample_steps):
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
    
    # sample a new sequence
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
            struct = random.choice(design.structures)
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

def classic_optimization(dg, design, objective_function=calculate_objective, weight=1, exit=1000, mode='sample', progress=False):
    '''
    Takes a Design object and does a classic optimization of this sequence.
    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_function: function which takes a design object and returns a score for evaluation
    :param weight: float specifying the weight of the difference part of the objective function
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

    score = objective_function(design, weight)
    # count for exit condition
    count = 0
    # remember how may mutations were done
    number_of_samples = 0
    
    # main optimization loop 
    while True:
        # count up the mutations
        number_of_samples += 1
        # sample a new sequence
        (mut_nos, sample_count) = _sample_sequence(dg, design, mode, 1)
        
        # write progress
        if progress:
            sys.stdout.write("\rMutate: {0:7.0f}/{1:5.0f} | Score: {2:7.4f} | NOS: {3:.5e}".format(number_of_samples, count, score, mut_nos) + " " * 20)
            sys.stdout.flush()
        
        this_score = objective_function(design, weight)
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
        sys.stdout.write("\r" + " " * 60 + "\r")
        sys.stdout.flush()
    
    # finally return the result
    return score, number_of_samples

def constraint_generation_optimization(dg, design, objective_function=calculate_objective, weight=1, exit=1000, mode='sample', num_neg_constraints=100, max_eos_diff=0, progress=False):
    '''
    Takes a Design object and does a constraint generation optimization of this sequence.
    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_function: function which takes a design object and returns a score for evaluation
    :param weight: float specifying the weight of the difference part of the objective function
    :param exit: Number of unsuccessful new sequences before exiting the optimization
    :param mode: String defining the sampling mode: sample, sample_global, sample_local
    :param num_neg_constraints: Maximal number of negative constraints to accumulate during the optimization process
    :param max_eos_diff: Maximal difference between eos of the negative and positive constraints
    :param progress: Whether or not to print the progress to the console
    :param return: Optimization score reached for the final sequence
    "param return: Number of samples neccessary to reach this result
    '''
    
    neg_constraints = collections.deque(maxlen=num_neg_constraints)
    
    # if the design has no sequence yet, sample one from scratch
    if not design.sequence:
        dg.sample()
        design.sequence = dg.get_sequence()
    else:
        dg.set_sequence(design.sequence)
    
    score = objective_function(design, weight);
    # count for exit condition
    count = 0
    # remember how may mutations were done
    number_of_samples = 0
    # sample steps to do
    sample_steps = 1
    # count for unsucessful constraint generation sequences
    cg_count = 0
    
    # main optimization loop
    while True:
        # constraint generation loop
        while True:
            # count up the mutations
            number_of_samples += 1
            # evaluate cg_count and make search space bigger if necessary
            if (cg_count > 10000):
                cg_count = 0
                sample_steps += 1
                dg.set_history_size(sample_steps+100)
            # increase cg_count
            cg_count += 1
            
            # sample a new sequence
            (mut_nos, sample_count) = _sample_sequence(dg, design, mode, sample_steps)
            
            # write progress
            if progress:
                sys.stdout.write("\rMutate: {0:7.0f}/{1:5.0f} | Steps: {2:3d} | EOS-Diff: {3:4.2f} | Score: {4:7.4f} | NOS: {5:.5e}".format(number_of_samples, count, sample_steps, max_eos_diff, score, mut_nos) + " " * 20)
                sys.stdout.flush()
            # boolean if it is perfect already
            perfect = True
            # evaluate the constraints
            for x in range(0, len(neg_constraints)):
                # test if the newly sampled sequence is compatible to the neg constraint, if not -> Perfect!
                if rd.sequence_structure_compatible(design.sequence, [neg_constraints[x]]):
                    if design.classtype == 'vrnaDesign':
                        neg_eos = RNA.energy_of_struct(design.sequence, neg_constraints[x])
                    elif design.classtype == 'nupackDesign':
                        neg_eos = nupack.energy([design.sequence], neg_constraints[x], material = 'rna', pseudo = True)
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
        # if we reached the mfe strcture, calculate a score for this solution and evaluate
        if design.mfe_structure in design.structures:
            # calculate objective
            this_score = objective_function(design, weight)
            # evaluate
            if (this_score < score):
                score = this_score
                # reset values
                sample_steps = 1
                cg_count = 0
                count = 0
            else:
                dg.revert_sequence(sample_count)
                design.sequence = dg.get_sequence()
        # else if current mfe is not in negative constraints, add to it
        else:
            if design.mfe_structure not in neg_constraints:
                neg_constraints.append(design.mfe_structure)
                #print('\n'+'\n'.join(neg_constraints))

            dg.revert_sequence(sample_count)
            design.sequence = dg.get_sequence()
        
        # exit condition
        if count > exit:
            break
        
    # clear the console
    if (progress):
        sys.stdout.write("\r" + " " * 60 + "\r")
        sys.stdout.flush()
    # finally return the result
    return score, number_of_samples

def _sample_connected_components(dg, amount):
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
    for k in range(solution_space_size - sample_size + 1, solution_space_size + 1):
        sum += float(solution_space_size) / float(k)
    return float(sum)

