#!/usr/bin/env python
'''
    test_State_Class.py: UNIT tests for State_Class.py
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
    
    def test_ensemble_defect(self):
        a = nupackDesign(['((((((((((....))))))))))'], 'GCCCCCCCCGGAAACGGGGGGGGC')
        ed = a.state['0'].ensemble_defect
        print('nupack: {0:4.7f}'.format(ed))
        self.assertEqual(ed, 0.0257900)
        
        b = vrnaDesign(['((((((((((....))))))))))'], 'GCCCCCCCCGGAAACGGGGGGGGC')
        ed = b.state['0'].ensemble_defect
        print('vrna ed: {0:4.7f}'.format(ed))
        self.assertEqual(round(ed, 7), 0.0592835)
        

if __name__ == '__main__':
    unittest.main()