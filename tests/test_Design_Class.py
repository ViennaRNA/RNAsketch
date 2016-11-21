#!/usr/bin/env python
'''
    test_Design_Class.py: UNIT tests for Design_Class.py
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from RNAsketch import *
import unittest

class TestDesignClass(unittest.TestCase):

    def test_init(self):
        a = vrnaDesign(['.'], 'A')
        self.assertEqual(a.sequence, 'A')
        self.assertEqual(a.structures, ['.'])
        with self.assertRaises(TypeError):
            b = vrnaDesign('.', 'A')
        with self.assertRaises(TypeError):
            b = vrnaDesign(['.'], 'Y')
        with self.assertRaises(TypeError):
            b = vrnaDesign(['%'], 'A')
        b = vrnaDesign({'eins': '..', 'zwei': '()'}, 'AU')
        self.assertEqual(b.state['eins'].structure, '..')
        self.assertEqual(b.state['zwei'].structure, '()')
    
    def test_init_empty_sequence(self):
        a = vrnaDesign(['...'])
        self.assertEqual(a.sequence, None)
        self.assertEqual(a.structures, ['...'])
        self.assertEqual(a.mfe_energy, {'0': None})
        self.assertEqual(a.mfe_structure, {'0': None})
        self.assertEqual(a.pf_energy, {'0': None})
        self.assertEqual(a.pf_structure, {'0': None})
        self.assertEqual(a.number_of_structures, 1)
        self.assertEqual(a.eos, {'0': None})
        self.assertEqual(a.pos, {'0': None})
        self.assertEqual(a.eos_reached_mfe, {'0': None})
        self.assertEqual(a.eos_diff_mfe, {'0': None})
        self.assertEqual(a.ensemble_defect, {'0': None})
        
        b = vrnaDesign({'first': '...'})
        self.assertEqual(b.structures, ['...'])
        self.assertEqual(b.state['first'].structure, '...')
    
    def test_calculations(self):
        a = vrnaDesign(['.'], 'A')
        self.assertEqual(a.mfe_energy, {'0': 0.0})
        self.assertEqual(a.mfe_structure, {'0': '.'})
        self.assertEqual(a.pf_energy, {'0': 0.0})
        self.assertEqual(a.pf_structure, {'0': '.'})
        self.assertEqual(a.number_of_structures, 1)
        self.assertEqual(a.length, 1)
        self.assertEqual(a.eos, {'0': 0.0})
        self.assertEqual(a.pos, {'0': 1.0})
        self.assertEqual(a.eos_reached_mfe, {'0': 1})
        self.assertEqual(a.eos_diff_mfe, {'0': 0.0})
        self.assertEqual(a.ensemble_defect, {'0': 0.0})
    
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
        assert a.mfe_energy is not {'0': None}
        assert a.eos is not {'0': None}
        assert a.mfe_energy != mfe_energy
        assert a.pf_energy != pf_energy

if __name__ == '__main__':
    unittest.main()
