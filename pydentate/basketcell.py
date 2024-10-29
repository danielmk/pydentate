# -*- coding: utf-8 -*-
"""
@author: DanielM
"""

from ouropy.genneuron import GenNeuron
import ouropy.parameters as params
import os


class BasketCell(GenNeuron):
    """Implements the BasketCell class with logic inherited from
    ouropy.GenNeuron, a generic neuron class"""
    name = "BasketCell"

    def __init__(self, name=None):
        # Make soma
        self.mk_soma(name='soma', diam=15, L=20)

        # Make dendrites
        self.mk_dendrite(4, dend_name='basal_1',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[50, 50, 50, 50], soma_loc=0)
        self.mk_dendrite(4, dend_name='basal_2',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[50, 50, 50, 50], soma_loc=0)
        self.mk_dendrite(4, dend_name='apical_1',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[75, 75, 75, 75], soma_loc=1.0)
        self.mk_dendrite(4, dend_name='apical_2',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[75, 75, 75, 75], soma_loc=1.0)

        dirname = os.path.dirname(__file__)
        filepath = os.path.join(dirname, 'basketcellparams.txt')
        parameters = params.read_parameters(filepath)

        self.insert_mechs(parameters)

class BasketCell_Ih(GenNeuron):
    """Implements the BasketCell class with logic inherited from
    ouropy.GenNeuron, a generic neuron class"""
    name = "BasketCell"

    def __init__(self, name=None):
        # Make soma
        self.mk_soma(name='soma', diam=15, L=20)

        # Make dendrites
        self.mk_dendrite(4, dend_name='basal_1',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[50, 50, 50, 50], soma_loc=0)
        self.mk_dendrite(4, dend_name='basal_2',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[50, 50, 50, 50], soma_loc=0)
        self.mk_dendrite(4, dend_name='apical_1',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[75, 75, 75, 75], soma_loc=1.0)
        self.mk_dendrite(4, dend_name='apical_2',
                         sec_names=['proxd', 'mid1', 'mid2', 'ddend'],
                         diam=[4, 3, 2, 1], L=[75, 75, 75, 75], soma_loc=1.0)

        dirname = os.path.dirname(__file__)
        filepath = os.path.join(dirname, 'basketcellparams_Ih.txt')
        parameters = params.read_parameters(filepath)

        self.insert_mechs(parameters)
