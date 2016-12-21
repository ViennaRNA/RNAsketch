#!/usr/bin/env python
'''
    test_State_Class.py: UNIT tests for State_Class.py
'''

__author__ = "Stefan Hammer"
__copyright__ = "Copyright 2016"
__maintainer__ = "Stefan Hammer"
__email__ = "s.hammer@univie.ac.at"

from RNAsketch import *
import unittest

class TestStateClass(unittest.TestCase):
    
    def test_init(self):
        a = vrnaState(None, '(...)')
        self.assertEqual(a._parent, None)
        self.assertEqual(a.structure, '(...)')
        self.assertEqual(a.temperature, 37.0)
        self.assertEqual(a.ligand, None)
        self.assertEqual(a.constraint, None)
        self.assertEqual(a.enforce_constraint, False)
        self.assertEqual(a.length, 5)
        self.assertEqual(a.cut_points, [])
        self.assertEqual(a.multifold, 0)
        
    def test_set_values(self):
        a = vrnaState(None, '(...)')
        a.temperature = 24.0
        self.assertEqual(a.temperature, 24.0)
        a.ligand = ['A&U', '(&)', -5.6]
        self.assertEqual(a.ligand, ['A&U', '(&)', -5.6])
        a.constraint = '.(.).'
        self.assertEqual(a.constraint, '.(.).')
        a.enforce_constraint = True
        self.assertEqual(a.enforce_constraint, True)
    
    def test_wrong_values(self):
        with self.assertRaises(ValueError):
            a = vrnaState(None, ['(...)'])
        a = vrnaState(None, '(...)')
        with self.assertRaises(AttributeError):
            a.structure = '.....'
        with self.assertRaises(ValueError):
            a.temperature = 't'
        with self.assertRaises(ValueError):
            a.ligand = ['A&U', '(&)']
        with self.assertRaises(ValueError):
            a.constraint = '.(.)'
        with self.assertRaises(ValueError):
            a.constraint = '.(.))'
    
    def test_cut_points(self):
        a = vrnaState(None, '..&..+..&..')
        self.assertEqual(a.cut_points, [3,6,9])
        a = vrnaState(None, '..&..+..&&..')
        self.assertEqual(a.cut_points, [3,6,9,10])
    
    def test_eos(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].mfe_energy, a.state['0'].eos)
        self.assertNotEqual(a.state['1'].mfe_energy, a.state['1'].eos)
    
    def test_pos(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(round(a.state['0'].pos, 2), 0.94)
        self.assertEqual(round(a.state['1'].pos, 2), 0.0)
    
    def test_eos_diff_mfe(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].eos_diff_mfe, 0)
        self.assertEqual(a.state['1'].eos_diff_mfe, a.state['1'].eos-a.state['1'].mfe_energy)
    
    def test_eos_reached_mfe(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].eos_reached_mfe, True)
        self.assertEqual(a.state['1'].eos_reached_mfe, False)
    
    def test_mfe_structure(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].mfe_structure, '((((....))))')
        self.assertEqual(a.state['1'].mfe_structure, '((((....))))')
    
    def test_mfe_energy(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].mfe_energy, -5.0)
        self.assertEqual(a.state['1'].mfe_energy, a.state['0'].mfe_energy)
    
    def test_pf_structure(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].pf_structure, '((((....))))')
        self.assertEqual(a.state['1'].pf_structure, a.state['0'].pf_structure)
    
    def test_pf_energy(self):
        a = vrnaDesign(['((((....))))','............'], 'CCGCAAAAGCGG')
        self.assertEqual(round(a.state['0'].pf_energy, 2), -5.04)
        self.assertEqual(a.state['1'].pf_energy, a.state['0'].pf_energy)
    
    def test_ensemble_defect(self):
        a = nupackDesign(['((((((((((....))))))))))'], 'GCCCCCCCCGGAAACGGGGGGGGC')
        ed = a.state['0'].ensemble_defect
        print('nupack: {0:4.7f}'.format(ed))
        self.assertEqual(ed, 0.0257900)
        
        b = vrnaDesign(['((((((((((....))))))))))'], 'GCCCCCCCCGGAAACGGGGGGGGC')
        ed = b.state['0'].ensemble_defect
        print('vrna ed: {0:4.7f}'.format(ed))
        self.assertEqual(round(ed, 7), 0.0592835)
    
    def test_temperature(self):
        a = vrnaDesign(['((((....))))','((((....))))'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].pf_energy, a.state['1'].pf_energy)
        a.state['1'].temperature = 24.0
        self.assertNotEqual(a.state['0'].pf_energy, a.state['1'].pf_energy)
    
    def test_ligand(self):
        a = vrnaDesign(['((((....))))','((((....))))'], 'CCGCAAAAGCGG')
        self.assertEqual(a.state['0'].pf_energy, a.state['1'].pf_energy)
        a.state['1'].ligand = ['CC&GG', '((&))', -5.6]
        self.assertNotEqual(a.state['0'].pf_energy, a.state['1'].pf_energy)
        self.assertEqual(round(a.state['0'].mfe_energy, 2), round(a.state['1'].mfe_energy+5.6, 2))
    
    def test_length(self):
        a = vrnaState(None, '(...)')
        self.assertEqual(a.length, 5)
        with self.assertRaises(AttributeError):
            a.length = 6
    
    def test_constraints(self):
        a = vrnaDesign(['((((((((((....))))))))))'], 'GCCCCCCCCGGAAACGGGGGGGGC')
        pf = a.state['0'].pf_energy
        a.state['0'].constraint = '........((....))........'
        cpf = a.state['0'].pf_energy
        a.state['0'].enforce_constraint = True
        ecpf = a.state['0'].pf_energy
        print('{0:4.7f}\t{1:4.7f}\t{2:4.7f}'.format(pf, cpf, ecpf))
        self.assertNotEqual(pf, cpf)
        self.assertNotEqual(cpf, ecpf)
        self.assertNotEqual(pf, ecpf)
        

if __name__ == '__main__':
    unittest.main()
