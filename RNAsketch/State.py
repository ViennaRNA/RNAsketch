#!/usr/bin/env python
'''
    State.py: Class as wrapper for ViennaRNA and Nupack
    functions to design an RNA molecule.
    This implements the various states of a riboswitch.
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"


import re
import math
import sys

vrna_available = True
nupack_available = True
pKiss_available = True
HotKnots_available = True

try:
    import RNA
    #RNA.read_parameter_file('/usr/share/ViennaRNA/rna_turner1999.par')
except ImportError as e:
    vrna_available = False
    sys.stderr.write("-" * 60 + "\nWARNING: " + str(e) + "!!!\n" + "-" * 60 + "\n")
    sys.stderr.flush()

try:
    import nupack
except ImportError as e:
    nupack_available = False
    sys.stderr.write("-" * 60 + "\nWARNING: " + str(e) + "!!!\n" + "-" * 60 + "\n")
    sys.stderr.flush()

try:
    import pKiss
except ImportError as e:
    pKiss_available = False
    sys.stderr.write("-" * 60 + "\nWARNING: " + str(e) + "!!!\n" + "-" * 60 + "\n")
    sys.stderr.flush()

try:
    import HotKnots
except ImportError as e:
    HotKnots_available = False
    sys.stderr.write("-" * 60 + "\nWARNING: " + str(e) + "!!!\n" + "-" * 60 + "\n")
    sys.stderr.flush()

if not (vrna_available or nupack_available or pKiss_available or HotKnots_available):
    raise ImportError("Neither ViennaRNA package, Nupack or pKiss found!")

class State(object):
    '''
    State object holds structure, ligand, constraints and calculates
    all values for this state.

    :param structure: Dot-bracket structure string
    :param parent: Parent Design object
    '''

    def __init__(self, parent, structure=None, temperature=37.0, ligand=None, constraint=None, enforce_constraint=False):
        if not isinstance(structure, basestring):
            raise ValueError('Not a structure string: ' + repr(structure))

        self._structure = structure
        self._parent = parent
        self.reset()
        self._temperature = temperature
        self._ligand = ligand
        self._constraint = constraint
        self._enforce_constraint = enforce_constraint
        self._length = None
        self._cut_points = None
        self._multifold = None

    def reset(self):
        self._eos = None
        self._pos = None
        self._eos_diff_mfe = None
        self._eos_reached_mfe = None
        self._mfe_structure = None
        self._mfe_energy = None
        self._pf_structure = None
        self._pf_energy = None
        self._ensemble_defect = None

    @property
    def structure(self):
        '''
        :return: Dot-bracket structure as constraint for this state
        '''
        return self._structure

    @property
    def temperature(self):
        '''
        :return: Temperature of this state
        '''
        return self._temperature
    @temperature.setter
    def temperature(self, t):
        if not isinstance(t, float):
            raise ValueError('Temperature must be a float pointing value')
        self.reset()
        self._temperature = t

    @property
    def ligand(self):
        '''
        :return: Ligand binding definitions for this state; List with sequence pattern, structure pattern and binding energy (soft constraint)
        '''
        return self._ligand
    @ligand.setter
    def ligand(self, lig):
        if (len(lig) != 3) or not isinstance(lig, list):
            raise ValueError('Ligand must be a list with three items: sequence pattern, structure pattern and binding energy')
        self.reset()
        self._ligand = lig

    @property
    def constraint(self):
        '''
        :return: Structural constraint of this state (hard constraint)
        '''
        return self._constraint
    @constraint.setter
    def constraint(self, constraint):
        if constraint:
            create_bp_table(constraint) #check for balanced brackets
            if (len(constraint) != self.length):
                raise ValueError('constraint and sequence must have equal length!')
        self.reset()
        self._constraint = constraint

    @property
    def enforce_constraint(self):
        '''
        :return: Boolean whether to enforce the hard constraint or not
        '''
        return self._enforce_constraint
    @enforce_constraint.setter
    def enforce_constraint(self, enforced):
        if enforced is not self._enforce_constraint:
            self.reset()
            self._enforce_constraint = enforced

    @property
    def length(self):
        '''
        :return: Length of the design sequence if available, otherwise of the stuctural input
        '''
        if not self._length:
            if self._parent and self._parent.sequence:
                self._length = len(self._parent.sequence)
            else:
                self._length = len(self._structure)

        return self._length

    @property
    def cut_points(self):
        '''
        :return: List containing the positions of the cut points (& or +) in the sequence
        '''
        if not self._cut_points and self._structure:
            self._cut_points = []
            iterator = re.finditer(re.compile('\&|\+'), self._structure)
            for match in iterator:
                self._cut_points.append(match.start()+1)
        return self._cut_points

    @property
    def multifold(self):
        '''
        :return: Number of cut points to specify if it is a cofold or multifold problem
        '''
        if not self._multifold:
            self._multifold = len(self.cut_points)
        return self._multifold

    @property
    def classtype(self):
        return None

    @property
    def eos(self):
        '''
        :return: Energy of structure given all the properties of this state (constraints, temperature,...)
        '''
        if not self._eos and self._parent.sequence and self._structure:
            self._eos = self._get_eos(self._parent.sequence, self._structure, self.temperature, self.ligand)
        return self._eos

    @property
    def pos(self):
        '''
        :return: Probability of structure given all the properties of this state (constraints, temperature,...)
        '''
        if not self._pos and self._parent.sequence:
            self._pos = math.exp((self.pf_energy-self.eos) / self._get_KT(self.temperature) )
        return self._pos

    @property
    def eos_diff_mfe(self):
        '''
        :return: Energy difference of structure to mfe given all the properties of this state (constraints, temperature,...)
        '''
        if not self._eos_diff_mfe and self._parent.sequence:
            self._eos_diff_mfe = self.eos - self.mfe_energy
        return self._eos_diff_mfe

    @property
    def eos_reached_mfe(self):
        '''
        :return: Boolean defining whether the eos reached the mfe energy
        '''
        if not self._eos_reached_mfe and self._parent.sequence:
            if (self.eos == self.mfe_energy):
                self._eos_reached_mfe = 1
            else:
                self._eos_reached_mfe = 0
        return self._eos_reached_mfe

    @property
    def mfe_energy(self):
        '''
        :return: MFE value given all the properties of this state (constraints, temperature,...)
        '''
        if not self._mfe_energy and self._parent.sequence:
            self._calculate_mfe_energy_structure()
        return self._mfe_energy

    @property
    def mfe_structure(self):
        '''
        :return: MFE structure given all the properties of this state (constraints, temperature,...)
        '''
        if not self._mfe_structure and self._parent.sequence:
            self._calculate_mfe_energy_structure()
        return self._mfe_structure

    def _calculate_mfe_energy_structure(self):
        (structure, energie) = self._get_fold(self._parent.sequence, self.temperature, self.ligand, self.constraint)
        self._mfe_energy = energie
        self._mfe_structure = structure

    @property
    def pf_energy(self):
        '''
        In case of a dimer cofold calculation with viennarna, the function returns the energy of true hybrid states only.

        :return: Partition function energy value given all the properties of this state (constraints, temperature,...)
        '''
        if not self._pf_energy and self._parent.sequence:
            self._calculate_pf_energy_structure()
        return self._pf_energy

    @property
    def pf_structure(self):
        '''
        :return: Partition function consensus structure given all the properties of this state (constraints, temperature,...)
        '''
        if not self._pf_structure and self._parent.sequence:
            self._calculate_pf_energy_structure()
        return self._pf_structure

    def _calculate_pf_energy_structure(self):
        (structure, energie) = self._get_pf_fold(self._parent.sequence, self.temperature, self.ligand, self.constraint)
        self._pf_energy = energie
        self._pf_structure = structure

    @property
    def ensemble_defect(self):
        '''
        :return: Ensemble defect value given all the properties of this state (constraints, temperature,...)
        '''
        if not self._ensemble_defect and self._parent.sequence and self._structure:
            if (len(self._parent.sequence) != len(self._structure)):
                raise ValueError('sequence and structure must have equal length to calculate the ensemble defect!')
            self._ensemble_defect = self._get_ensemble_defect(self._parent.sequence, self._structure, self.temperature, self.ligand)
        return self._ensemble_defect

    def _get_KT(self, temperature):
        # KT = (betaScale*((temperature+K0)*GASCONST))/1000.0; /* in Kcal */
        return ((temperature + 273.15)*1.98717)/1000.0;

    def _get_eos(self, sequence, structure, temperature, ligand):
        raise NotImplementedError

    def _get_fold(self, sequence, temperature, ligand, constraint):
        raise NotImplementedError

    def _get_pf_fold(self, sequence, temperature, ligand, constraint):
        raise NotImplementedError

    def _get_ensemble_defect(self, sequence, structure, temperature, ligand):
        raise NotImplementedError

if vrna_available:
    class vrnaState(State):
        @property
        def classtype(self):
            return 'vrna'

        def _change_cuts(self, input):
            return re.sub('[+]', '&', input)

        def _get_fold_compound(self, sequence, temperature, ligand=None, constraint=None, options=RNA.OPTION_PF):
            md = RNA.md()
            md.temperature = temperature
            md.dangles = 2
            fc = RNA.fold_compound(self._change_cuts(sequence), md, options)
            if ligand:
                fc.sc_add_hi_motif(ligand[0], ligand[1], ligand[2])
            if constraint:
                if self._enforce_constraint:
                    fc.hc_add_from_db(remove_cuts(constraint), RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)
                else:
                    fc.hc_add_from_db(remove_cuts(constraint))
            return fc

        def _get_eos(self, sequence, structure, temperature, ligand=None):
            if self.multifold > 1:
                raise NotImplementedError
            fc = self._get_fold_compound(sequence, temperature, ligand, options=RNA.OPTION_MFE | RNA.OPTION_EVAL_ONLY)
            return fc.eval_structure(remove_cuts(structure))

        def _get_fold(self, sequence, temperature, ligand=None, constraint=None):
            fc = self._get_fold_compound(sequence, temperature, ligand, constraint, options=RNA.OPTION_MFE)
            if self.multifold == 0:
                (structure, energie) = fc.mfe()
            if self.multifold == 1:
                (structure, energie) = fc.mfe_dimer()
                structure = add_cuts(structure, self.cut_points)
            if self.multifold > 1:
                raise NotImplementedError
            return (structure, energie)

        def _get_pf_fold(self, sequence, temperature, ligand=None, constraint=None):
            fc = self._get_fold_compound(sequence, temperature, ligand, constraint)
            if self.multifold == 0:
                (structure, energie) = fc.pf()
            if self.multifold == 1:
                # returns (string structure, float *FA, float *FB, float *FcAB, float *FAB)
                (structure, _, _, energie, _) = fc.pf_dimer()
                structure = add_cuts(structure, self.cut_points)
            elif self.multifold > 1:
                raise NotImplementedError
            return (structure, energie)

        def _get_ensemble_defect(self, sequence, structure, temperature, ligand=None):
            #TODO finish implementation with pairing matrix. need to ask ronny how it is done now or use old interface
            fc = self._get_fold_compound(sequence, temperature, ligand)
            if self.multifold == 0:
                _,_ = fc.pf()
            if self.multifold == 1:
                _,_ = fc.pf_dimer()
                structure = add_cuts(structure, self.cut_points)
            elif self.multifold > 1:
                raise NotImplementedError
            # get base pairing probability matrix
            bpm = fc.bpp()
            # get base pair table
            bpt = create_bp_table(structure)
            # delta(i,j, s) = 1 if (i,j) in s, 0 otherwise
            # d(s1, s2) = sum{i,j in s1}(1-delta{i,j}(s2)) + sum{i,j not in s1}(delta{i,j}(s2))
            # ensemble defect = sum{i,j in structure}(1-P(i,j)) + sum{i,j not in structure}(P(i,j))
            result = 0
            for i in range(0, len(structure)):
                for j in range(i, len(structure)):
                    if (bpt[i] == j):
                        #print('{0:1.0f}/{1:1.0f}: {2:10.10f}'.format(i, j, bpm[i+1][j+1]))
                        result += 1-bpm[i+1][j+1]
                    else:
                        result += bpm[i+1][j+1]
            return 2 * result

if nupack_available:
    class nupackState(State):
        @property
        def classtype(self):
            return 'nupack'

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

        def _get_ensemble_defect(self, sequence, structure, temperature, ligand=None):
            return nupack.defect([self._change_cuts(sequence)], structure, material = 'rna', pseudo = True, T = temperature)

if pKiss_available:
    class pKissState(State):
        @property
        def classtype(self):
            return 'pKiss'

        def _change_cuts(self, input):
            if re.match(r'[&+]', input):
                raise IOError('pKiss cannot handle concatenated RNAs')
            return input

        def _get_eos(self, sequence, structure, temperature, ligand=None):
            #TODO pKiss cannot handle ligands
            # pkiss eval returns multiple answers to the eval question, take lowest energy!
            shape, energy, structure = pKiss.eval(self._change_cuts(sequence), self._change_cuts(structure), temperature = temperature)
            return energy

        def _get_fold(self, sequence, temperature, ligand=None, constraint=None):
            return pKiss.mfe(self._change_cuts(sequence), temperature = temperature)

        def _get_pf_fold(self, sequence, temperature, ligand=None, constraint=None):
            #TODO return mfe as a bad approximation it cannot calculate the partition function
            return pKiss.mfe(self._change_cuts(sequence), temperature = temperature)

        def _get_ensemble_defect(self, sequence, structure, temperature, ligand=None):
            raise NotImplementedError

if HotKnots_available:
    class hotknotsState(State):
        @property
        def classtype(self):
            return 'hotknots'

        def _change_cuts(self, input):
            if re.match(r'[&+]', input):
                raise IOError('Hotknots cannot handle concatenated RNAs')
            return input

        def _get_eos(self, sequence, structure, temperature, ligand=None):
            #TODO Hotknots cannot handle ligands and temperature
            structure, energy = HotKnots.eval(self._change_cuts(sequence), self._change_cuts(structure))
            return energy

        def _get_fold(self, sequence, temperature, ligand=None, constraint=None):
            #TODO Hotknots cannot handle ligands and temperature
            return HotKnots.mfe(self._change_cuts(sequence))

        def _get_pf_fold(self, sequence, temperature, ligand=None, constraint=None):
            #TODO return mfe as a bad approximation it cannot calculate the partition function
            return HotKnots.mfe(self._change_cuts(sequence))

        def _get_ensemble_defect(self, sequence, structure, temperature, ligand=None):
            raise NotImplementedError

def remove_cuts(input):
    '''
    Takes a string and removes all cut characters (+&) from this string

    :param input: string containing cut point characters
    :return: string without cut point characters
    '''
    return re.sub('[+&]', '', input)

def add_cuts(input, cut_points):
    '''
    Takes a string and add back the cut point characters at the positions stored in the cut points list

    :param input: string without cut point characters
    :return: string containing cut point characters
    '''
    result = input
    for cut in cut_points:
        result = result[:cut-1] + '&' + result[cut-1:]
    return result

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
        raise ValueError('Unbalanced brackets: too few closing brackets')
    return bpt
