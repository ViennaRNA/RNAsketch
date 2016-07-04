#!/usr/bin/env python
'''
    testPyDesign.py: UNIT tests for Design_Class.py
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from PyDesign import *
import unittest

class TestDesignClass(unittest.TestCase):

    def test_init(self):
        a = vrnaDesign(['.'], 'A')
        self.assertEqual(a.sequence, 'A')
        self.assertEqual(a.structures, ['.'])
        with self.assertRaises(TypeError):
            b = vrnaDesign('.', 'A')
            b = vrnaDesign(['.'], 'Y')
            b = vrnaDesign(['%'], 'A')
            
    def test_init_empty_sequence(self):
        a = vrnaDesign(['...'])
        self.assertEqual(a.sequence, None)
        self.assertEqual(a.structures, ['...'])
        self.assertEqual(a.mfe_energy, [None])
        self.assertEqual(a.mfe_structure, [None])
        self.assertEqual(a.pf_energy, [None])
        self.assertEqual(a.pf_structure, [None])
        self.assertEqual(a.number_of_structures, 1)
        self.assertEqual(a.eos, [None])
        self.assertEqual(a.pos, [None])
        self.assertEqual(a.eos_reached_mfe, [None])
        self.assertEqual(a.eos_diff_mfe, [None])
      

    def test_calculations(self):
        a = vrnaDesign(['.'], 'A')
        self.assertEqual(a.mfe_energy, [0.0])
        self.assertEqual(a.mfe_structure, ['.'])
        self.assertEqual(a.pf_energy, [0.0])
        self.assertEqual(a.pf_structure, ['.'])
        self.assertEqual(a.number_of_structures, 1)
        self.assertEqual(a.length, 1)
        self.assertEqual(a.eos, [0.0])
        self.assertEqual(a.pos, [1.0])
        self.assertEqual(a.eos_reached_mfe, [1])
        self.assertEqual(a.eos_diff_mfe, [0])
    
    def test_print(self):
        a = vrnaDesign(['((((....))))','..((....))..'], 'AAGGACGUCCUU')
        print '\n'
        print a.write_out(1000)
        print '\n'
        print a.write_csv_header()
        print a.write_csv()
        print '\n'
    
    def test_reassign(self):
        a = vrnaDesign(['((((....))))'], 'AAAAGGGGUUUU')
        mfe_energy = a.mfe_energy
        pf_energy = a.pf_energy
        a.sequence = 'GGGGAAAACCCC'
        assert a.mfe_energy is not None
        assert a.eos is not None
        assert a.mfe_energy != mfe_energy
        assert a.pf_energy != pf_energy