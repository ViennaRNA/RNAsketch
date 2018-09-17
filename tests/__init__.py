#!/usr/bin/env python
'''
    test_RNAsketch.py: UNIT tests for RNAsketch.py
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from RNAsketch import *
from test_State import TestStateClass
from test_Design import TestDesignClass
import tempfile
import unittest

class TestRNAsketch(unittest.TestCase):

    def test_read_inp_file(self):
        fp = tempfile.NamedTemporaryFile(delete=True)
        fp.write(b'.....\n  C  \n;')
        fp.seek(0)
        structures, constraint, sequence = read_inp_file(fp.name)
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
        structures, constraint, sequence = read_input(data)
        self.assertEqual(structures, ['.[].','({})'])
        self.assertEqual(constraint, 'NNAC')
        self.assertEqual(sequence, 'AGCC')

    def test_read_input_additions(self):
        data = '.[].;;5.,;5\n({})\t\nAGCC\t\t\nNNAC'
        structures, constraint, sequence, additions = read_input_additions(data)
        self.assertEqual(structures, ['.[].','({})'])
        self.assertEqual(additions, ['5.,;5',''])
        self.assertEqual(constraint, 'NNAC')
        self.assertEqual(sequence, 'AGCC')

    def test_get_graph_properties(self):
        pass #this is not necessary to test

    def test_calculate_objective(self):
        a = vrnaDesign(['((((....))))'], 'GGGGAAAACCCC')
        score = calculate_objective(a)
        self.assertEqual(score, 0.107421875)

    def test_adaptive_walk_optimization(self):
        pass

    def test_constraint_generation_optimization(self):
        pass

    def test_sample_sequence(self):
        pass

    def _sample_connected_components(self):
        pass

    def test_sample_count_unique_solutions(self):
        self.assertEqual(sample_count_unique_solutions(6, 6), 14.7)



if __name__ == '__main__':

    unittest.main()
