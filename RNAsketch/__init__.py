#!/usr/bin/env python
'''
    RNAsketch.py: A small module containing all the important helpers and
    functions for all major RNAdesign operations.
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

import sys
import random
import collections
import numpy as np
import math

import RNAblueprint as rbp
from Design import *

'''
Global variable:
'''

def read_inp_file(filename):
    '''
    Reads a file in .inp format and returns all neccessary information

    :param filename: Filename of the file to read
    :return: structures - List of structures in dot-bracket notation
    :return: constraint - Sequence constraint
    :return: sequence - Start sequence
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
    :return: structures - List of structures in dot-bracket notation
    :return: constraint - Sequence constraint
    :return: sequence - Start sequence
    '''
    return read_input_additions(content)[:3]

def read_input_additions(content):
    '''
    Reads some input and returns all neccessary information in the right container.
    Input is a string, lines separated by linebreaks. Content might be structures,
    a sequence constraint and a start sequence. Additional information for the structural
    states can be provided with any separator ;,: or whitespaces after the structure.

    :param filename: Filename of the file to read
    :return: structures - List of structures in dot-bracket notation
    :return: constraint - Sequence constraint
    :return: sequence - Start sequence
    :return: additions - List of additional information after the structures
    '''
    structures = []
    constraint = ''
    sequence = ''
    additions = []

    lines = content.split("\n")
    for line in lines:
        # if line begins with a semicolon ; stop parsing
        if re.match(re.compile("^\;"), line, flags=0):
            break
        # strip additional information after the structure/sequence string
        m = re.match(re.compile("^([^\s\;\,\:]+)[\s\;\,\:]*(.*)$"), line, flags=0)
        if m:
            line, addition = m.groups()
        if re.match(re.compile("^[\(\)\.\{\}\[\]\<\>\+\&]+$"), line, flags=0):
            structures.append(line.rstrip('\n'))
            additions.append(addition.rstrip('\n'))
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

    return structures, constraint, sequence, additions

def get_graph_properties(dg):
    '''
    Takes a RNAdesign DependencyGraph Object and constructs a dicionary with all the
    calculated properties.

    :param dg: RNAdesign DependencyGraph object
    :return: properties - Dictionary containing all the graph properties
    '''
    properties = {}
    articulation_ratios = []
    max_articulations = 0
    max_component_vertices = 0

    properties['num_cc'] = dg.number_of_connected_components()
    properties['nos'] = dg.number_of_sequences()

    for cc in range(0, properties['num_cc']):
        cv = len(dg.component_vertices(cc))
        sv = len(dg.articulation_vertices(cc))
        articulation_ratios.append(float(sv)/float(cv))
        if (max_articulations < sv):
            max_articulations = sv
        if (max_component_vertices < cv):
            max_component_vertices = cv

    properties['max_articulations'] = max_articulations
    properties['max_component_vertices'] = max_component_vertices
    properties['max_articulation_ratio'] = max(articulation_ratios)
    properties['mean_articulation_ratio'] = sum(articulation_ratios) / len(articulation_ratios)
    properties['max_dimensions'] = dg.max_number_of_dimensions()

    return properties

def calculate_objective(design, weight=0.5):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    (eos(1)+eos(2)+eos(3) - 3 * gibbs) / number_of_structures +
    weight * (eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2) * 2 / (number_of_structures * (number_of_structures-1))

    :param design: Design object containing the sequence and structures
    :type design: Object of type Design
    :param weight: To wheight the influence of the eos diffences
    :type weight: float
    :return: score calculated by the objective function
    '''
    return calculate_objective_1(design) + weight * calculate_objective_2(design)

def calculate_objective_1(design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    (eos(1)+eos(2)+eos(3) - 3 * gibbs) / number_of_structures

    :param design: Design object containing the sequence and structures
    :type design: Object of type Design
    :return: score calculated by the objective function
    '''
    return (sum(design.eos.values()) - sum(design.pf_energy.values())) / design.number_of_structures

def calculate_objective_2(design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    |eos(1)-eos(2)| + |eos(1)-eos(3)| + |eos(2)-eos(3))| * 2 / (number_of_structures * (number_of_structures-1))

    :param design: Design object containing the sequence and structures
    :return: score calculated by the objective function
    '''
    objective_difference_part = 0
    eos = design.eos.values()
    for i, eos1 in enumerate(eos):
        for eos2 in eos[i+1:]:
            objective_difference_part += math.fabs(eos1 - eos2)

    if design.number_of_structures == 1:
        return objective_difference_part
    else:
        return objective_difference_part * 2 / (design.number_of_structures * (design.number_of_structures-1))

def calculate_objective_2_squared(design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    (eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2) * 2 / (number_of_structures * (number_of_structures-1))

    :param design: Design object containing the sequence and structures
    :return: score calculated by the objective function
    '''
    objective_difference_part = 0
    eos = design.eos.values()
    for i, eos1 in enumerate(eos):
        for eos2 in eos[i+1:]:
            objective_difference_part += (eos1 - eos2) ** 2

    if design.number_of_structures == 1:
        return objective_difference_part
    else:
        return objective_difference_part * 2 / (design.number_of_structures * (design.number_of_structures-1))


def sample_sequence(dg, design, mode, sample_steps=1, avoid_motifs=None, white_positions=None):
    '''
    This function samples a sequence with the given mode from the dependency graph object
    and writes it into the design object

    :param dg: RNAdesign dependency graph object
    :param design: design object
    :param mode: mode how to sample, this is a string
    :param sample_steps: count how many times to do the sample operation
    :param avoid_motifs: list of regex pattern specifiying sequence motifs to avoid
    :param white_positions: list of [start, end] positions in the sequence where the avoid_motifs pattern should be ignored
    :return: mut_nos is the solution space we drew from
    :return: sample_count is how many times we sampled a solution from the dependency graph object (important for revert later)
    '''
    if avoid_motifs is None:
        avoid_motifs=[]
    if white_positions is None:
        white_positions=[]
    # remember the solution space we drew from
    mut_nos = 1
    dg.set_history_size(sample_steps + 100)
    while True:
        # count how many samples we did to be able to revert this later
        sample_count = 0
        # sample a new sequence
        # set the chosen mode for this run
        chosen_mode = mode
        # if random choice is requested pick something new
        if mode == "random":
            modes = ['sample','sample_clocal','sample_plocal']
            chosen_mode = random.choice(modes)

        if sample_steps == 0:
            sample_steps = random.randrange(1, dg.number_of_connected_components())

        if chosen_mode == 'sample':
            mut_nos = dg.sample()
            sample_count += 1
        elif chosen_mode == 'sample_clocal':
            for c in _sample_connected_components(dg, sample_steps):
                mut_nos *= dg.sample_clocal(c)
                sample_count += 1
        elif chosen_mode == 'sample_plocal':
            # TODO this local sampling is unfair this way and mut_nos is not calculated correctly!
            for _ in range(0, sample_steps):
                mut_nos *= dg.sample_plocal()
                sample_count += 1
        elif consen_mode == 'sample_unpaired':
            # sample only unpaired positions
            for _ in range(0, sample_steps):
                mut_nos *= dg.sample_clocal(1,1)
                sample_count += 1
        else:
            raise ValueError("Wrong mode argument: " + mode + "\n")

        # check if motifs to avoid are present, if so sample a new sequence, else return
        motiv_present = False
        seq = dg.get_sequence()
        # search for all motivs and ignore all hits inside any white_position
        for m in avoid_motifs:
            for f in re.finditer(re.compile(r'(?=('+m+'))'), seq, flags=0):
                motiv_present = True
                for white in white_positions:
                    if ((white[0] <= f.start() <= white[1]) and (white[0] <= f.start()+len(f.group(1))-1 <= white[1])):
                        motiv_present = False
                        break
                if motiv_present:
                    break
            if motiv_present:
                break

        # if the motivs are not present, exit while and return
        if not motiv_present:
            break
        else:
            # revert to previous sequence
            dg.revert_sequence(sample_count)

    # assign sequence to design and return values
    design.sequence = dg.get_sequence()
    return (mut_nos, sample_count)

def adaptive_walk_optimization(dg, design, objective_function=calculate_objective, stop=1000, mode='sample', avoid_motifs=None, white_positions=None, progress=False):
    '''
    Takes a Design object and does a adaptive walk optimization of this sequence.

    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_functions: array of functions which takes a design object and returns a score for evaluation
    :param stop: Number of unsuccessful new sequences before stoping the optimization
    :param mode: String defining the sampling mode: sample, sample_clocal, sample_plocal
    :param avoid_motifs: list of regex pattern specifiying sequence motifs to avoid
    :param white_positions: list of [start, end] positions in the sequence where the avoid_motifs pattern should be ignored
    :param progress: Whether or not to print the progress to the console
    :return: Optimization score reached for the final sequence
    :return: Number of samples neccessary to reach this result
    '''
    if avoid_motifs is None:
        avoid_motifs=[]
    if white_positions is None:
        white_positions=[]
    # if the design has no sequence yet, sample one from scratch
    if not design.sequence:
        sample_sequence(dg, design, 'sample', avoid_motifs=avoid_motifs, white_positions=white_positions)
    else:
        dg.set_sequence(design.sequence)

    score = objective_function(design)
    # count for stop condition
    count = 0
    # remember how may mutations were done
    number_of_samples = 0

    # main optimization loop
    while stop:
        # count up the mutations
        number_of_samples += 1
        # sample a new sequence
        (mut_nos, sample_count) = sample_sequence(dg, design, mode, avoid_motifs=avoid_motifs, white_positions=white_positions)

        # write progress
        if progress:
            sys.stderr.write("\rMutate: {0:7.0f}/{1:5.0f} | Score: {2:5.2f} | NOS: {3:.5e} | Mode: {4:s}".format(number_of_samples, count, score, mut_nos, mode) + " " * 20)
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
            if count > stop:
                break

    # clear the console
    if (progress):
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()
    # finally return the result
    return score, number_of_samples

def simulated_annealing_optimization(dg, design, objective_function=calculate_objective, temperature_gradient=None, cooling_step=50, mode='sample', avoid_motifs=None, white_positions=None, progress=False):
    '''
    Takes a Design object and does a simulated annealing optimization of this sequence.

    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_functions: array of functions which takes a design object and returns a score for evaluation
    :param temperature_gradient: Iterable containing the temperatures in descend order
    :param cooling_steps: Use the current temperature that many times before cooling down
    :param mode: String defining the sampling mode: sample, sample_clocal, sample_plocal
    :param avoid_motifs: list of regex pattern specifiying sequence motifs to avoid
    :param white_positions: list of [start, end] positions in the sequence where the avoid_motifs pattern should be ignored
    :param progress: Whether or not to print the progress to the console
    :return: Optimization score reached for the final sequence
    :return: Number of samples neccessary to reach this result
    '''
    if avoid_motifs is None:
        avoid_motifs=[]
    if white_positions is None:
        white_positions=[]
    if temperature_gradient is None:
        temperature_gradient=np.concatenate([np.arange(1,0,-0.0002),[1e-15]*100])
    # generate iterator (can call next() on it)
    temp_iter = iter(temperature_gradient)
    temperature = temp_iter.next()

    # if the design has no sequence yet, sample one from scratch
    if not design.sequence:
        sample_sequence(dg, design, 'sample', avoid_motifs=avoid_motifs, white_positions=white_positions)
    else:
        dg.set_sequence(design.sequence)

    score = objective_function(design)
    # remember how may mutations were done
    number_of_samples = 0
    # remember how often we used the temperature already
    number_of_same_temp = 0

    # main optimization loop
    while True:
        # check if we have to cool down
        number_of_same_temp += 1
        if number_of_same_temp > cooling_step:
            number_of_same_temp = 0
            try:
                temperature = temp_iter.next()
            except StopIteration:
                # end of temperature scale reached... stop optimization
                break

        # count up the mutations
        number_of_samples += 1
        # sample a new sequence
        (mut_nos, sample_count) = sample_sequence(dg, design, mode, avoid_motifs=avoid_motifs, white_positions=white_positions)

        # write progress
        if progress:
            sys.stderr.write("\rMutate: {0:7.0f}/{1:5.0f} | Score: {2:5.2f} | NOS: {3:.5e} | Mode: {4:s} | Temp: {5:5.8f}".format(number_of_samples, number_of_same_temp, score, mut_nos, mode, temperature) + " " * 20)
            sys.stderr.flush()

        this_score = objective_function(design)
        # evaluate probability
        rand = random.uniform(0, 1)
        if (this_score-score) < 0:
            prob = 1
        else:
            prob = math.exp(-1*(this_score-score)/temperature)
        # compare and make decision
        if (rand <= prob):
            score = this_score
            # go to next temperature step
            number_of_same_temp = cooling_step+1
        else:
            dg.revert_sequence(sample_count)
            design.sequence = dg.get_sequence()

    # clear the console
    if (progress):
        sys.stderr.write("\r" + " " * 60 + "\r")
        sys.stderr.flush()
    # finally return the result
    return score, number_of_samples

def constraint_generation_optimization(dg, design, objective_function=calculate_objective, stop=1000, mode='sample', num_neg_constraints=100, max_eos_diff=0, avoid_motifs=None, white_positions=None, progress=False):
    '''
    Takes a Design object and does a constraint generation optimization of this sequence.

    :param dg: RNAdesign DependencyGraph object
    :param design: Design object containing the sequence and structures
    :param objective_functions: array of functions which takes a design object and returns a score for evaluation
    :param stop: Number of unsuccessful new sequences before stoping the optimization
    :param mode: String defining the sampling mode: sample, sample_clocal, sample_plocal
    :param num_neg_constraints: Maximal number of negative constraints to accumulate during the optimization process
    :param max_eos_diff: Maximal difference between eos of the negative and positive constraints
    :param avoid_motifs: list of regex pattern specifiying sequence motifs to avoid
    :param white_positions: list of [start, end] positions in the sequence where the avoid_motifs pattern should be ignored
    :param progress: Whether or not to print the progress to the console
    :return: Optimization score reached for the final sequence
    :return: Number of samples neccessary to reach this result
    '''
    if avoid_motifs is None:
        avoid_motifs=[]
    if white_positions is None:
        white_positions=[]
    dg.set_history_size(100)
    neg_constraints = collections.deque(maxlen=num_neg_constraints)

    # if the design has no sequence yet, sample one from scratch
    if not design.sequence:
        sample_sequence(dg, design, 'sample', avoid_motifs=avoid_motifs, white_positions=white_positions)
    else:
        dg.set_sequence(design.sequence)

    score = objective_function(design)
    # count for stop condition
    count = 0
    # remember how may mutations were done
    number_of_samples = 0

    # main optimization loop
    while stop:
        # constraint generation loop
        while True:
            # count up the mutations
            number_of_samples += 1
            # sample a new sequence
            (mut_nos, sample_count) = sample_sequence(dg, design, mode, avoid_motifs=avoid_motifs, white_positions=white_positions)

            # write progress
            if progress:
                sys.stderr.write("\rMutate: {0:7.0f}/{1:5.0f} | EOS-Diff: {2:4.2f} | Score: {3:5.2f} | NOS: {4:.5e} | Mode: {5:s}".format(number_of_samples, count, max_eos_diff, score, mut_nos, mode))
                sys.stderr.flush()
            # boolean if it is perfect already
            perfect = True
            # evaluate the constraints
            for negc in reversed(neg_constraints):
                # test if the newly sampled sequence is compatible to the neg constraint, if not -> Perfect!
                if rbp.sequence_structure_compatible(design.sequence, [negc]):
                    if design.classtype == 'vrna':
                        neg_eos = RNA.energy_of_struct(design.sequence, negc)
                    elif design.classtype == 'nupack':
                        neg_eos = nupack.energy([design.sequence], negc, material = 'rna', pseudo = True)
                    else:
                        raise ValueError('Could not figure out the classtype of the Design object.')
                    # test if the newly sampled sequence eos for pos constraints is lower than
                    # the eos for all negative constraints, if not -> Perfect!
                    for eos in design.eos.values():
                        if float(neg_eos) - float(eos) < max_eos_diff:
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
        for mfe_str in design.mfe_structure.values():
            if mfe_str not in design.structures:
                if mfe_str not in neg_constraints:
                    neg_constraints.append(mfe_str)
                    #print('\n'+'\n'.join(neg_constraints))

        # stop condition
        if count > stop:
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
    :return: list of connected component IDs which can be used for example for: dg.sample_clocal(ID)
    '''
    result = []
    noslist = {}

    for c in range(0, dg.number_of_connected_components()):
        #if dg.number_of_sequences(c) != 1:
        noslist[c] = dg.number_of_sequences(c)

    if amount > len(noslist):
        amount = len(noslist)

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
