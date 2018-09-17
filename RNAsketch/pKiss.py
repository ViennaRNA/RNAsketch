#!/usr/bin/env python
from __future__ import print_function

'''
    pKiss.py: Wrapper for the pKiss software, a tool for folding RNA secondary structures, also ones with two limited classes of pseudoknots.
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2018"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

import math
import subprocess as sub
import os
import tempfile
import re

def _setup_args(**kargs):
    '''
    Takes the arguments from the various calls and generates a string array for the command-line call.

    :param kargs: List of the function arguments
    :return: List of strings with the command line arguments
    :return: String with command line input
    '''

    mode = kargs['mode']
    args =  ['--mode', kargs['mode']]

    if mode != 'abstract':
        args += ['--temperature', kargs['temperature'],
            '--allowLP', int(kargs['lonelyPairs'])]


    if mode in ["subopt", "local", "shapes", "cast"]:
        if kargs['relativeDeviation']:
            # deviation in percent from mfe
            args += ['--relativeDeviation', kargs['deviation']]
        else:
            # deviation in kcal/mol from mfe
            args += ['--absoluteDeviation', kargs['deviation']]

    if mode in ["shapes", "probs", "cast", "eval", "abstract"]:
        args += ['--shapeLevel', kargs['shapeLevel']] # 1-5 most abstract is 5

    if 'window' in kargs and kargs['window']:
        if mode in ["mfe", "subopt", "enforce", "local", "shapes", "probs"]:
            args += ['--windowSize', kargs['windowSize'], # >0 do the calculation for every window with this size
            '--windowIncrement', kargs['windowIncrement']] # next window with this stepsize

    if mode in ["mfe", "subopt", "enforce", "local", "shapes", "probs", "cast"]:
        args += ['--strategy', kargs['pkStrategy'], # pseudoknot strategy, A B C D P, default: A
        '--maxKnotSize', kargs['pkMaxKnot']] # max pseudoknot size, default: Length

    if mode in ["mfe", "subopt", "enforce", "local", "shapes", "probs", "cast", "eval"]:
            args += ['--minHairpinLength', kargs['pkMinHairpin'] # length of initial hairpin, less better but slower
            ]
            #'--Hpenalty', kargs['Hpenalty'], # openin H-type pseudoknot energy, default 9kcal/mol
            #'--Kpenalty', kargs['Kpenalty']] # openin K-type pseudoknot energy, default 12kcal/mol

    if mode == "probs":
        args += ['--lowProbFilter', kargs['discardStates']] # 0-<1filter results with low probabilities, default 0
        #'--outputLowProbFilter' # 0-<1 only discard on output

    # command line arguments must be strings
    args = list(map(str, args))
    if mode == 'eval':
        cmd_input = '\n'.join(['>1', kargs['sequence'], kargs['structure']])
    elif mode == 'abstract':
        cmd_input = kargs['structure']
    else:
        # stdin is just the sequence
        cmd_input = kargs['sequence']

    return args, cmd_input

def _filter_errors(error):
    rvalues = []

    if error:
        pattern = re.compile(r'Redundant argument', flags=0)
        error_lines = error.decode().rstrip('\n').split('\n')
        for l in error_lines:
            if not re.match(pattern, l):
                rvalues.append(l)
    return '\n'.join(rvalues)

def _call_with_pipe(args, cmd_input):
    '''
    Performs the command line call.

    :param args: List of strings containing the program arguments
    :param cmd_input: String containing the program stdin input
    :return: String with output of the program
    '''
    #print(" ".join(['# pKiss'] + args))
    p=sub.Popen(['pKiss'] + args, stdin = sub.PIPE, stdout = sub.PIPE, stderr = sub.PIPE)
    output, error = p.communicate(cmd_input)

    e = _filter_errors(error)
    if e:
        raise IOError('pKiss returned an error:\n' + e)

    return output.decode()

def _call_with_file(args, cmd_input):
    '''
    Command line call with a tmp file as input.

    :param args: List of strings containing the program arguments
    :param cmd_input: String containing the program stdin input
    :return: String with output of the program
    '''
    inputfile = tempfile.NamedTemporaryFile(delete=False, suffix=".fa")

    inputfile.write(cmd_input.encode())
    inputfile.close()

    #print(" ".join(['# pKiss'] + args))
    p = sub.Popen(['pKiss'] + args + [inputfile.name], stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.STDOUT)
    output, error = p.communicate(None)

    e = _filter_errors(error)
    if e:
        raise IOError('pKiss returned an error:\n' + e)

    os.remove(inputfile.name)
    return output.decode()

def _parse_output(output, multistate = False, window = False):
    '''
    Regex that parses the output of the programs.

    :param output: Output string of the pKiss program
    :param multistate: Boolean stating if it is a multistate output
    :param window: Boolean stating if it is window mode output
    :return: Either a tuple of structure/shape and energy values, or in window mode a dictionary with the start position as key and the tuple as value.
    '''
    #print(output)
    pattern = re.compile(r'(?P<start>\d+)\s+(?P<sequence>[AUGC]+)\s+(?P<end>\d+)\n^\s*(?P<states>[-\d\s\[\]\{\}\(\)\<\>\.\_]+\n)', flags = re.M)
    state_pattern = re.compile(r'(?P<energy>-?\d+\.?\d*)\s*(?P<struct>[\[\]\{\}\(\)\<\>\.]+)(\s+(?P<prob>\d+\.\d+)?\s+(?P<shape>[\[\]\{\}\(\)\<\>\.\_]+))?', flags = 0)

    values = {}
    for m in re.finditer(pattern, output):
        values[int(m.group('start'))] = []
        #print(m.group('states'))
        for s in re.finditer(state_pattern, m.group('states')):
            if s.group('shape'):
                if s.group('prob'):
                    values[int(m.group('start'))].append( (s.group('shape'), float(s.group('prob')), s.group('struct'), float(s.group('energy')) ) )
                else:
                    values[int(m.group('start'))].append( (s.group('shape'), float(s.group('energy')), s.group('struct') ) )
            else:
                # struct and energy must always be present
                values[int(m.group('start'))].append( (s.group('struct'), float(s.group('energy')) ) )
        # flatten states list if only one
        if not multistate:
            values[int(m.group('start'))] = values[int(m.group('start'))][0]
    # flatten dict if no window mode
    if not values:
        raise IOError('Something went wrong calling pKiss:\n' + output)
    if not window:
        return values[1]
    else:
        return values

def mfe(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None):
    '''
    Computes the single energetically most stable secondary structure for the
    given RNA sequence. This structure might contain a pseudoknot of type H
    (simple canonical recursive pseudoknot) or type K (simple canonical
    recursive kissing hairpin), but need not to. Co-optimal results will be
    suppressed, i.e. should different prediction have the same best energy
    value, just an arbitrary one out of them will be reported.

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :return: Tuple of string and float containing the mfe structure and energy
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'mfe', sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=False, window=False)

def mfe_window(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None, windowSize = 5, windowIncrement = 1):
    '''
    Window mode of mfe.
    Computes the single energetically most stable secondary structure for the
    given RNA sequence. This structure might contain a pseudoknot of type H
    (simple canonical recursive pseudoknot) or type K (simple canonical
    recursive kissing hairpin), but need not to. Co-optimal results will be
    suppressed, i.e. should different prediction have the same best energy
    value, just an arbitrary one out of them will be reported.

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :param windowSize: Int size of the shifting window (default: 5)
    :param windowIncrement: Int jump size of the shifting window (default: 1)
    :return: Dictionary with start position as key and as key the tuple of string and float containing the mfe structure and energy
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'mfe', window = True, sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot, windowSize = windowSize, windowIncrement = windowIncrement)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=False, window=True)

def subopt(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None, relativeDeviation = True, deviation=20.0):
    '''
    Often, the biological relevant structure is hidden among suboptimal
    predictions. In "subopt mode", you can also inspect all suboptimal solutions
    up to a given threshold (see parameters --absoluteDeviation and
    --relativeDeviation). Due to semantic ambiguity of the underlying
    "microstate" grammar, sometimes identical predictions will show up. As
    Vienna-Dot-Bracket strings they seem to be the same, but according to base
    dangling they differ and thus might even have slightly different energies.
    See 1] for details.

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :param relativeDeviation: Boolean stating whether the energy deviation should is given in percent of the mfe energy (default: True)
    :param deviation: Float energy deviation above the mfe energy either in kcal/mol or in percent of the mfe energy (default: 20%)
    :return: List containing the states in tuples of string and float representing the structure and energy
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'subopt', sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot, relativeDeviation = relativeDeviation, deviation = deviation)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=True, window=False)

def subopt_window(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None, relativeDeviation = True, deviation=20.0, windowSize = 5, windowIncrement = 1):
    '''
    Often, the biological relevant structure is hidden among suboptimal
    predictions. In "subopt mode", you can also inspect all suboptimal solutions
    up to a given threshold (see parameters --absoluteDeviation and
    --relativeDeviation). Due to semantic ambiguity of the underlying
    "microstate" grammar, sometimes identical predictions will show up. As
    Vienna-Dot-Bracket strings they seem to be the same, but according to base
    dangling they differ and thus might even have slightly different energies.
    See 1] for details.

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :param relativeDeviation: Boolean stating whether the energy deviation should is given in percent of the mfe energy (default: True)
    :param deviation: Float energy deviation above the mfe energy either in kcal/mol or in percent of the mfe energy (default: 20%)
    :param windowSize: Int size of the shifting window (default: 5)
    :param windowIncrement: Int jump size of the shifting window (default: 1)
    :return: Dictionary with start position as key and as key the list containing the states in tuples of string and float representing the structure and energy
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'subopt', window = True, sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot, relativeDeviation = relativeDeviation, deviation = deviation, windowSize = windowSize, windowIncrement = windowIncrement)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=True, window=True)

def enforce():
    '''
    Energetically best pseudoknots might be deeply buried under suboptimal
    solutions. Use "enforce" mode to enforce a structure prediction for each of
    the for classes: "nested structure" (as "RNAfold" would compute, i.e.
    without pseudoknots), "H-type pseudoknot", "K-type pseudoknot" and "H- and
    K-type pseudoknot". Useful if you want to compute the tendency of folding a
    pseudoknot or not, like in 2].
    '''
    raise NotImplementedError

def local():
    '''
    Computes energetically best and suboptimal local pseudoknots. Local means,
    leading and trailing bases can be omitted and every prediction is a
    pseudoknot.
    '''
    raise NotImplementedError

def shapes(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None, relativeDeviation = True, deviation=20.0, shapeLevel = 2):
    '''
    Output of "subopt" mode is crowded by many very similar answers, which make
    it hard to focus to the "important" changes. The abstract shape concept 3]
    groups similar answers together and reports only the best answer within such
    a group. Due to abstraction, suboptimal analyses can be done more thorough,
    by ignoring boring differences. (see parameter --shapeLevel)

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :param relativeDeviation: Boolean stating whether the energy deviation should is given in percent of the mfe energy (default: True)
    :param deviation: Float energy deviation above the mfe energy either in kcal/mol or in percent of the mfe energy (default: 20%)
    :param shapeLevel: Int abstraction level of the shape representation (default: 2)
    :return: List containing the states in tuples of string, float and string representing the shape, energy and structure.
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'shapes', sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot, relativeDeviation = relativeDeviation, deviation = deviation, shapeLevel = shapeLevel)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=True, window=False)

def shapes_window(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None, relativeDeviation = True, deviation=20.0, shapeLevel = 2, windowSize = 5, windowIncrement = 1):
    '''
    Window version of shapes.
    Output of "subopt" mode is crowded by many very similar answers, which make
    it hard to focus to the "important" changes. The abstract shape concept 3]
    groups similar answers together and reports only the best answer within such
    a group. Due to abstraction, suboptimal analyses can be done more thorough,
    by ignoring boring differences. (see parameter --shapeLevel)

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :param relativeDeviation: Boolean stating whether the energy deviation should is given in percent of the mfe energy (default: True)
    :param deviation: Float energy deviation above the mfe energy either in kcal/mol or in percent of the mfe energy (default: 20%)
    :param shapeLevel: Int abstraction level of the shape representation (default: 2)
    :param windowSize: Int size of the shifting window (default: 5)
    :param windowIncrement: Int jump size of the shifting window (default: 1)
    :return: Dictionary with start position as key and as key the list containing the states in tuples of string, float and string representing the shape, energy and structure
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'shapes', window = True, sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot, relativeDeviation = relativeDeviation, deviation = deviation, shapeLevel = shapeLevel, windowSize = windowSize, windowIncrement = windowIncrement)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=True, window=True)

def probs(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None, shapeLevel = 2, discardStates = 0.0):
    '''
    Structure probabilities are strictly correlated to their energy values.
    Grouped together into shape classes, their probabilities add up. Often a
    shape class with many members of worse energy becomes more probable than the
    shape containing the mfe structure but not much more members. See 4] for
    details on shape probabilities.

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :param discardStates: Float [0,1) probility treshold under which states should be discarded during the calculation. (default: 0.0)
    :param shapeLevel: Int abstraction level of the shape representation (default: 2)
    :return: List containing the states in tuples of string, float, string and float representing the shape, probability, structure and energy
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'probs', sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot, shapeLevel = shapeLevel, discardStates = discardStates)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=True, window=False)

def probs_window(sequence, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, pkMaxKnot = None, shapeLevel = 2, discardStates = .1, windowSize = 5, windowIncrement = 1):
    '''
    Window version of probs.
    Structure probabilities are strictly correlated to their energy values.
    Grouped together into shape classes, their probabilities add up. Often a
    shape class with many members of worse energy becomes more probable than the
    shape containing the mfe structure but not much more members. See 4] for
    details on shape probabilities.

    :param sequence: String containing the RNA sequence [AUGC]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param pkMaxKnot: Int maximal PK size (default: Length of sequence)
    :param discardStates: Float [0,1) probility treshold under which states should be discarded during the calculation. (default: 0.0)
    :param shapeLevel: Int abstraction level of the shape representation (default: 2)
    :param windowSize: Int size of the shifting window (default: 5)
    :param windowIncrement: Int jump size of the shifting window (default: 1)
    :return: Dictionary with start position as key and as key the list containing the states in tuples of string, float, string and float representing the shape, probability, structure and energy
    '''
    if pkMaxKnot is None:
        pkMaxKnot = len(sequence)

    args, cmd_input = \
    _setup_args(mode = 'probs', window = True, sequence = sequence, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, pkMaxKnot = pkMaxKnot, shapeLevel = shapeLevel, discardStates = discardStates, windowSize = windowSize, windowIncrement = windowIncrement)

    output = _call_with_pipe(args, cmd_input)

    return _parse_output(output, multistate=True, window=True)

def cast():
    '''
    This mode is the RNAcast approache, see 8]. For a family of RNA sequences,
    this method independently enumerates the near-optimal abstract shape space,
    and predicts as the consensus an abstract shape common to all sequences. For
    each sequence, it delivers the thermodynamically best structure which has
    this common shape. Input is a multiple fasta file, which should contain at
    least two sequences. Output is sorted by "score" of common shapes, i.e.
    summed free energy of all sequences. R is the rank (= list position) of the
    shape in individual sequence analysis.
    '''
    raise NotImplementedError

def eval(sequence, structure, temperature = 37.0, lonelyPairs = True, pkStrategy = 'A', pkMinHairpin = 2, shapeLevel = 2):
    '''
    Evaluates the free energy of an RNA molecule in fixed secondary structure,
    similar to RNAeval from the Vienna group. Multiple answers stem from
    semantic ambiguity of the underlying grammar. It might happen, that your
    given structure is not a structure for the sequence. Maybe your settings are
    too restrictive, e.g. not allowing lonely base-pairs (--allowLP), too long
    hairpin stems for pseudoknots (--minHairpinLength) or the given pseudoknot
    is more complex than those of pKiss. If you input a (multiple) FASTA file,
    pKiss assumes that exactly first half of the contents of each entry is RNA
    sequence, second half is the according structure. Whitespaces are ignored.

    :param sequence: String containing the RNA sequence [AUGC]
    :param structure: String containing the dot-bracket structure [.(){}[]<>]
    :param temperature: Float temperature in degree celsius (default: 37.0)
    :param lonelyPairs: Boolean stating whether lonely pairs are allowed (default: True)
    :param pkStrategy: String stating the pseudoknot evaluation strategy [A,B,C,D,P] (default: A)
    :param pkMinHairpin: Integer length of initial PK hairpin, less is better but slower. (default: 2)
    :param discardStates: Float [0,1) probility treshold under which states should be discarded during the calculation. (default: 0.0)
    :param shapeLevel: Int abstraction level of the shape representation (default: 2)
    :return: Tuple of string, float and string representing the shape, energy and structure
    '''
    if len(structure) != len(sequence):
        raise IOError('Sequence and Structure have unequal length!')

    args, cmd_input = \
    _setup_args(mode = 'eval', sequence = sequence, structure = structure, temperature = temperature, lonelyPairs = lonelyPairs, pkStrategy = pkStrategy, pkMinHairpin = pkMinHairpin, shapeLevel = shapeLevel)

    output = _call_with_file(args, cmd_input)

    return _parse_output(output, multistate=False, window=False)

def abstract(structure, shapeLevel = 2):
    '''
    Converts a Vienna-Dot-Bracket representation of a secondary structure into a
    shape string.

    :param structure: String containing the dot-bracket structure [.(){}[]<>]
    :param shapeLevel: Int abstraction level of the shape representation (default: 2)
    :return: String representing the shape
    '''

    args, cmd_input = \
    _setup_args(mode = 'abstract', structure = structure, shapeLevel = shapeLevel)

    output = _call_with_pipe(args, cmd_input)

    return output.rstrip('\n')

# Test for binary precense
try:
    output = _call_with_pipe(['--help'], '')
except:
    raise ImportError('pKiss not found. Please install pKiss!')
