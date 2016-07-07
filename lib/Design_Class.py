#!/usr/bin/env python
'''
    Design_Class.py: Class as wrapper for ViennaRNA and Nupack
    functions to design an RNA molecule
    This implements the core class which holds all the states of a riboswitch
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"


import re
from State_Class import *

class Design(object):
    '''
    Design contains all neccessary information to generate a RNA device
    or to return a solution. It is used as a container between the different
    functions.
    '''
    def __init__(self, structures, sequence=''):
        '''
        Construct a new State object.
        
        :param structures:
        :param sequence:
        '''
        self._number_of_structures = None
        self.state = {}
        self._structures = []
        
        if isinstance(structures, list):
            for key , struct in enumerate(structures):
                self._parseStructures(key, struct)
        elif isinstance(structures, dict):
            for key, struct in structures.items():
                self._parseStructures(key, struct)
        else:
            raise TypeError('Structures must be a list or a dict hoding the state name and the structure')
        
        self.sequence = sequence
        
    def _parseStructures(self, key, struct):
        create_bp_table(struct) #check for balanced brackets
        if not (isinstance(struct, basestring) and re.match(re.compile("[\(\)\.\+\&]"), struct)):
            raise TypeError('Structure be a string in dot-bracket notation')
        self.state[str(key)] = self._newState(struct)
        self._structures.append(struct)
    
    def _newState(self, struct):
        raise NotImplementedError
    
    @property
    def structures(self):
        return self._structures
    
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, s):
        if isinstance(s, basestring) and re.match(re.compile("[AUGC\+\&]"), s):
            self._sequence = s
            for state in self.state.values():
                state.reset()
        elif s == '':
            self._sequence = None
        else:
            raise TypeError('Sequence must be a string containing a IUPAC RNA sequence')
    
    @property
    def number_of_structures(self):
        if not self._number_of_structures:
            self._number_of_structures = len(self.structures)
        return self._number_of_structures
    
    def write_out(self, score=0):
        '''
        Generates a nice human readable version of all values of this design
        :param score: optimization score for this design
        :return: string containing a nicely formatted version of all design values
        '''
        result = '{0:}\t {1:5.2f}'.format(self.sequence, score)
        for k in self.state:
            state = self.state[k]
            result += '\n{0:}'.format(k)
            result += '\n{0:}\t{1:9.4f}\t{2:+9.4f}\t{3:9.4f}'.format(state.structure, state.eos, state.eos-state.mfe_energy, state.pos)
            result += '\n{0:}\t{1:9.4f}'.format(state.mfe_structure, state.mfe_energy)
            result += '\n{0:}\t{1:9.4f}'.format(state.pf_structure, state.pf_energy)
        return result
    
    def write_csv(self, separator=';'):
        '''
        Generates a csv version of all values of this design separated by the given separator
        :param separator: separator for the values
        :return: string containing all values of this design separated by the given separator
        '''
        result = separator.join(map(str, ['\"' + self.sequence + '\"', self.length, self.number_of_structures ]))
        for state in self.state.values():
            result = separator.join(map(str, [result,
            state.mfe_energy,
            state.mfe_structure,
            state.pf_energy,
            state.pf_structure,
            state.eos,
            state.eos_diff_mfe, 
            state.eos_reached_mfe,
            state.pos]))
        return result
    
    def write_csv_header(self, separator=';'):
        '''
        Generates a csv header for all values of this design separated by the given separator
        :param separator: separator for the values
        :return: string containing a csv header for this design separated by the given separator
        '''
        result = separator.join(['sequence', 'seq_length', 'number_of_structures'])
        strings = ['mfe_energy_', 'mfe_structure_', 'pf_energy_', 'pf_structure_', 'eos_', 'diff_eos_mfe_', 'mfe_reached_', 'prob_']
        for state in self.state:
            for s in strings:
                result += separator + s + state
        return result
    
    @property
    def eos(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].eos
        return result
    
    @property
    def pos(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].pos
        return result
    
    @property
    def eos_diff_mfe(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].eos_diff_mfe
        return result
    
    @property
    def eos_reached_mfe(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].eos_reached_mfe
        return result
    
    @property
    def mfe_structure(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].mfe_structure
        return result
    
    @property
    def mfe_energy(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].mfe_energy
        return result
    @property
    def pf_structure(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].pf_structure
        return result
    @property
    def pf_energy(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].pf_energy
        return result
    
    @property
    def ensemble_defect(self):
        result = {}
        for s in self.state:
            result[s] = self.state[s].ensemble_defect
        return result
    
    @property
    def length(self):
        return len(self.sequence)

class vrnaDesign(Design):
    def _newState(self, struct):
        return vrnaState(struct, self)

class nupackDesign(Design):
    def _newState(self, struct):
        return nupackState(struct, self)
