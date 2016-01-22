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

'''
Global variable:
'''
# KT = (betaScale*((temperature+K0)*GASCONST))/1000.0; /* in Kcal */
KT = ((37+273.15)*1.98717)/1000.0;

class Design:
    '''
    Design contains all neccessary information to generate a RNA device
    or to return a solution. It is used as a container between the different
    functions.
    '''
    def __init__(self, sequence, structures):
        '''
        Construct a new Design object.

        :param sequence:
        :param structures:
        '''
        self.sequence = sequence
        self.structures = structures
        self.score = None
        self.number_of_mutations = None
        self._mfe_energy
        self._mfe_structure
    
    @property
    def sequence(self):
        return self._sequence
    @sequence.setter
    def sequence(self, s):
        self._sequence = s
        self._reset_values()
    
    @property
    def structures(self):
        return self._structures
    @structures.setter
    def structures(self, s):
        self._structures = s
        self._reset_values()
        
    def _reset_values(self):
        self._eos = None
        self._mfe_structure = None
        self._mfe_energy = None
    
    @property
    def eos(self):
        if not self._eos:
            for struct in self.structures:
                self._eos.append(RNA.energy_of_struct(self.sequence, struct))
        return self._eos
    
    @property
    def mfe_energy(self):
        if not self._mfe_energy:
            (self._mfe_structure, self._mfe_energy) = RNA.fold(self.sequence)
        return self._mfe_energy
    
    @property
    def mfe_structure(self):
        if not self._mfe_structure:
            (self._mfe_structure, self._mfe_energy) = RNA.fold(self.sequence)
        return self._mfe_structure
    
    @property
    def partition_function(self):
        if not self._mfe_structure:
            self._part_funct = RNA.pf_fold(self.sequence)[1]
        return self._partition_function
    
    def number_of_structures(self):
        return len(self.structures)
    
        
    self.probs.append( math.exp((self.part_funct-this_eos) / kT ) )
    
    def write_out(self):
        #first clean up last line
        sys.stdout.write("\r" + " " * 60 + "\r")
        sys.stdout.flush()
        print(self.sequence + '\t{0:9.4f}'.format(self.score))
        for i, struct in enumerate(self.structures):
            print(struct + '\t{0:9.4f}\t{1:+9.4f}\t{2:9.4f}'.format(self.eos[i], self.eos[i]-self.mfe_energy, self.probs[i]))
        print(self.mfe_struct + '\t{0:9.4f}'.format(self.mfe_energy))


# objective function: eos(1)+eos(2)+eos(3) - 3 * gibbs + 1 * ((eos(1)-eos(2))^2 + (eos(1)-eos(3))^2 + (eos(2)-eos(3))^2)
def calculate_objective(Design):
    '''
    Calculates the objective function given a Design object containing the designed sequence and input structures.

    :param Design: Design object containing the sequence and structures
    '''
    eos = []
    for struct in Design.structures:
        eos.append(RNA.energy_of_struct(Design.sequence, struct))
    
    gibbs = RNA.pf_fold(Design.sequence)
    
    objective_difference_part = 0
    for i, value in enumerate(eos):
        for j in eos[i+1:]:
            objective_difference_part += math.fabs(value - j)
    
    number_structures = len(eos)
    return sum(eos) - number_structures * gibbs[1] + \
        1 * (objective_difference_part / (math.factorial(number_structures) / (math.factorial(number_structures-2)*2)))
