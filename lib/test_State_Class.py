#!/usr/bin/env python
'''
    testPyDesign.py: UNIT tests for State_Class.py
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__version__ = "0.1"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from PyDesign import *
import unittest

class TestStateClass(unittest.TestCase):
    
    
    def test_cut_points(self):
        a = vrnaDesign(['..&..+..&..'], 'AA&AA+AA')
        self.assertEqual(a.state['0'].cut_points, [3,6,9])
        a = vrnaDesign(['..&..+..&..'], 'AA&AA+AA')