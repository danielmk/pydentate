# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 14:44:15 2017

@author: DanielM
"""

import unittest

class TestNeuron(unittest.TestCase):
    """
    Imports a test neuron as specified by ouropy.tests.testneuron 
    and the parameters in ouropy/tests/testneuronparams.txt
    The tests check that all parameters were set as defined.
    """
    def setUp(self):
        """Create a generic neuron"""
        import ouropy.tests.testneuron
        self.testneuron = ouropy.tests.testneuron.TestNeuron()
        
    def test_soma_general(self):
        pass

    def test_soma_geometry(self):
        self.assertEqual(self.testneuron.soma.L, 16.8)
        self.assertEqual(self.testneuron.soma.diam, 20)
    
    def test_dendrites_general(self):
        self.assertEqual(len(self.testneuron.dendrites),2)
        self.assertEqual(len(self.testneuron.dendrites[0]), 4)
        self.assertEqual(len(self.testneuron.dendrites[1]), 4)
        
    def test_dendrites_geometry(self):
        length_dends = [50.0, 150.0, 150.0, 150.0]
        diam_dends = [5,4,3,2]
        for dend in self.testneuron.dendrites:
            for seg_idx, seg in enumerate(dend):
                self.assertEqual(seg.L, length_dends[seg_idx])
                self.assertEqual(seg.diam, diam_dends[seg_idx])
    
    def test_soma_biophysics(self):
        pass
        
class TestGenNeuronMethods(unittest.TestCase):
    """
    Imports the genneuron module and tests the methods of the GenNeuron class.
    """
    def setUp(self):
        from ouropy.genneuron import GenNeuron
        self.neuron = GenNeuron()
        
    def test_mk_soma(self):
        self.neuron.mk_soma(15,20)
        self.assertEqual(self.neuron.soma.diam, 15)
        self.assertEqual(self.neuron.soma.L, 20)
        self.neuron.mk_soma(20,25)
        self.assertEqual(self.neuron.soma.diam, 20)
        self.assertEqual(self.neuron.soma.L, 25)
        
    def test_new(self):
        print(self.neuron.)
        


    
        

if __name__ == '__main__':
    unittest.main()