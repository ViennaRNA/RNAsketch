#!/usr/bin/env python
'''
    testPyDesign.py: UNIT tests for PyDesign.py
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from PyDesign import *
import tempfile
import unittest

class TestPyDesign(unittest.TestCase):

    def test_init(self):
        a = Design(['.'], 'A')
        self.assertEqual(a.sequence, 'A')
        self.assertEqual(a.structures, ['.'])
        with self.assertRaises(TypeError):
            b = Design('.', 'A')
            b = Design(['.'], 'Y')
            b = Design(['%'], 'A')
    
    def test_init_empty_sequence(self):
        a = Design(['...'])
        self.assertEqual(a.sequence, None)
        self.assertEqual(a.structures, ['...'])
        self.assertEqual(a.mfe_energy, None)
        self.assertEqual(a.mfe_structure, None)
        self.assertEqual(a.pf_energy, None)
        self.assertEqual(a.pf_structure, None)
        self.assertEqual(a.number_of_structures, 1)
        self.assertEqual(a.eos, None)
        self.assertEqual(a.pos, None)
        self.assertEqual(a.eos_reached_mfe, None)
        self.assertEqual(a.eos_diff_mfe, None)
      

    def test_calculations(self):
        a = Design(['.'], 'A')
        self.assertEqual(a.mfe_energy, 0.0)
        self.assertEqual(a.mfe_structure, '.')
        self.assertEqual(a.pf_energy, 0.0)
        self.assertEqual(a.pf_structure, '.')
        self.assertEqual(a.number_of_structures, 1)
        self.assertEqual(a.length, 1)
        self.assertEqual(a.eos, [0.0])
        self.assertEqual(a.pos, [1.0])
        self.assertEqual(a.eos_reached_mfe, [1])
        self.assertEqual(a.eos_diff_mfe, [0])
    
    def test_print(self):
        a = Design(['((((....))))','..((....))..'], 'AAGGACGUCCUU')
        print '\n'
        print a.write_out(1000)
        print '\n'
        print a.write_csv_header()
        print a.write_csv()
        print '\n'
    
    def test_reassign(self):
        a = Design(['((((....))))'], 'AAAAGGGGUUUU')
        mfe_energy = a.mfe_energy
        pf_energy = a.pf_energy
        a.sequence = 'GGGGAAAACCCC'
        assert a.mfe_energy is not None
        assert a.eos is not None
        assert a.mfe_energy != mfe_energy
        assert a.pf_energy != pf_energy
    
    def test_read_inp_file(self):
        fp = tempfile.NamedTemporaryFile(delete=True)
        fp.write(b'.....\n  C  \n;')
        fp.seek(0)
        (structures, constraint, sequence) = read_inp_file(fp.name)
        self.assertEqual(structures, ['.....'])
        self.assertEqual(constraint, 'NNCNN')
        self.assertEqual(sequence, '')
    
    def test_wrong_inp_file(self):
        fp = tempfile.NamedTemporaryFile(delete=True)
        fp.write(b'....\n...\nNNY\n;')
        fp.seek(0)
        with self.assertRaises(ValueError):
            read_inp_file(fp.name)
    
    def test_read_input(self):
        data = '.[].\n({})\nAGCC\nNNAC'
        (structures, constraint, sequence) = read_input(data)
        self.assertEqual(structures, ['.[].','({})'])
        self.assertEqual(constraint, 'NNAC')
        self.assertEqual(sequence, 'AGCC')
    
    def test_get_graph_properties(self):
        print('TODO')
    
    def test_calculate_objective(self):
        a = Design(['((((....))))'], 'GGGGAAAACCCC')
        score = calculate_objective(a)
        self.assertEqual(score, 0.107421875)
    
    def test_classic_optimization(self):
        print('TODO')
        

if __name__ == '__main__':
    unittest.main()
