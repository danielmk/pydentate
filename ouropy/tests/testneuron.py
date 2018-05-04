# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 14:56:08 2017

@author: DanielM
"""

import ouropy.genneuron
import ouropy.parameters as params

class TestNeuron(ouropy.genneuron.GenNeuron):
    
    def __init__(self):
        self.mk_soma(name='soma', diam=20, L=16.8)

        # Make dendrites
        self.mk_dendrite(4, dend_name='dend_1',
                         sec_names=['gcld', 'proxd', 'midd', 'dd'],
                         diam=[5,4,3,2], L=[50.0, 150.0, 150.0, 150.0], soma_loc = 0)
        self.mk_dendrite(4, dend_name='dend_2',
                         sec_names=['gcld', 'proxd', 'midd', 'dd'],
                         diam=[5,4,3,2], L=[50.0, 150.0, 150.0, 150.0], soma_loc = 1)

        parameters = params.read_parameters('testneuronparams.txt')

        self.insert_mechs(parameters)
