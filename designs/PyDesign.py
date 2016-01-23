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
import RNAdesign as rd
import RNA
import math
import re

'''
Global variable:
'''
# KT = (betaScale*((temperature+K0)*GASCONST))/1000.0; /* in Kcal */
KT = ((37+273.15)*1.98717)/1000.0;

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
        self._sequence = None
        if sequence != '':
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
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, s):
        if isinstance(s, basestring) and re.match(re.compile("[AUGC]"), s):
            if len(s) != self.length:
                raise TypeError('Sequence must have the same length as the structural constraints')
            self._reset_values()
            self._sequence = s
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
                self._eos.append(RNA.energy_of_struct(self.sequence, struct))
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
            (self._mfe_structure, self._mfe_energy) = RNA.fold(self.sequence)
        return self._mfe_energy
    
    @property
    def mfe_structure(self):
        if not self._mfe_structure and self._sequence:
            (self._mfe_structure, self._mfe_energy) = RNA.fold(self.sequence)
        return self._mfe_structure
    
    @property
    def pf_energy(self):
        if not self._pf_energy and self._sequence:
            (self._pf_structure, self._pf_energy) = RNA.pf_fold(self.sequence)
        return self._pf_energy
    
    @property
    def pf_structure(self):
        if not self._pf_structure and self._sequence:
            (self._pf_structure, self._pf_energy) = RNA.pf_fold(self.sequence)
        return self._pf_structure
    
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
    
    def write_out(self, score=''):
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
            '\"'+ self.mfe_structure + '\"', self.pf_energy, '\"' + self.pf_structure + '\"']))
        result += separator + separator.join(map(str, self.eos))
        result += separator + separator.join(map(str, self.eos_diff_mfe))
        result += separator + separator.join(map(str, self.eos_reached_mfe))
        result += separator + separator.join(map(str, self.pos))
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


def calculate_objective(Design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.
    objective function (3 seqs):    eos(1)+eos(2)+eos(3) - 3 * gibbs + 
                                    weight * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2) / (3!/(3-2)!*2)
    :param Design: Design object containing the sequence and structures
    '''
    weight = 1
    
    objective_difference_part = 0
    eos = Design.eos
    for i, eos1 in enumerate(eos):
        for eos2 in eos[i+1:]:
            objective_difference_part += math.fabs(eos1 - eos2)
    
    combination_count = 1
    if (Design.number_of_structures != 1):
        combination_count = math.factorial(Design.number_of_structures) / (math.factorial(Design.number_of_structures-2)*2)
    
    return sum(Design.eos) - Design.number_of_structures * Design.pf_energy + weight * (objective_difference_part / combination_count)

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
        if re.match(re.compile("[\(\)\.]"), line, flags=0):
            structures.append(line.rstrip('\n'))
        elif re.match(re.compile("[\ ACGTUWSMKRYBDHVN]"), line, flags=0):
            line = line.replace(" ", "N")
            if re.match(re.compile("[ACGTU]"), line, flags=0) and sequence == '':
                sequence = line.rstrip('\n')
            elif constraint == '':
                constraint = line.rstrip('\n')
            else:
                raise ValueError('Too many constraints or start sequences')
    
    checklength = len(structures[0])
    if constraint != '' and len(constraint) != checklength:
        raise ValueError('Structures and the sequence constraint must have the same length!')
    elif sequence != '' and len(sequence) != checklength:
        raise ValueError('Structures and the start sequence must have the same length!')
    
    for s in structures:
        if len(s) != checklength:
            raise ValueError('Structures must all have the same length!')
    
    return structures, constraint, sequence