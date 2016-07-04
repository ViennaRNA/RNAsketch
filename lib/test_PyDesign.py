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
        a = vrnaDesign(['((((....))))'], 'GGGGAAAACCCC')
        score = calculate_objective(a)
        self.assertEqual(score, 0.107421875)
    
    def test_classic_optimization(self):
        print('TODO')
        

if __name__ == '__main__':
    unittest.main()
