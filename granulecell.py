# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:36:11 2017

@author: DanielM
"""
from ouropy.genneuron import GenNeuron
import ouropy.parameters as params


class GranuleCell(GenNeuron):
    """Implements the GranuleCell class with logic inherited from
    ouropy.GenNeuron, a generic neuron class"""
    name = 'GranuleCell'

    def __init__(self, name=None):
        # Make soma
        self.mk_soma(name='soma', diam=16.8, L=16.8)

        # Make dendrites
        self.mk_dendrite(4, dend_name='dend_1',
                         sec_names=['gcld', 'proxd', 'midd', 'dd'],
                         diam=[3.0, 3.0, 3.0, 3.0],
                         L=[50.0, 150.0, 150.0, 150.0], soma_loc=1.0)
        self.mk_dendrite(4, dend_name='dend_2',
                         sec_names=['gcld', 'proxd', 'midd', 'dd'],
                         diam=[3.0, 3.0, 3.0, 3.0],
                         L=[50.0, 150.0, 150.0, 150.0], soma_loc=1.0)

        parameters = params.read_parameters('granulecellparams.txt')

        self.insert_mechs(parameters)
