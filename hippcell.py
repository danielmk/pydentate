# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 10:42:49 2017

@author: DanielM
"""

from ouropy.genneuron import GenNeuron
import ouropy.parameters as params


class HippCell(GenNeuron):
    """Implements the HippCell class with logic inherited from
    ouropy.GenNeuron, a generic neuron class"""
    name = "HippCell"

    def __init__(self, name=None):
        self.mk_soma(name='soma', diam=10, L=20)

        self.mk_dendrite(3, dend_name='short_1',
                         sec_names=['proxd', 'midd', 'distd'],
                         diam=[3, 2, 1], L=[50, 50, 50], soma_loc=0.0)
        self.mk_dendrite(3, dend_name='short_2',
                         sec_names=['proxd', 'midd', 'distd'],
                         diam=[3, 2, 1], L=[50, 50, 50], soma_loc=0.0)
        self.mk_dendrite(3, dend_name='long_1',
                         sec_names=['proxd', 'midd', 'distd'],
                         diam=[3, 2, 1], L=[75, 75, 75], soma_loc=1.0)
        self.mk_dendrite(3, dend_name='long_2',
                         sec_names=['proxd', 'midd', 'distd'],
                         diam=[3, 2, 1], L=[75, 75, 75], soma_loc=1.0)

        parameters = params.read_parameters('hippcellparams.txt')

        self.parameters = parameters

        self.insert_mechs(parameters)
