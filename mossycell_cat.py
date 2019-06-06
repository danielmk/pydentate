# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:16:40 2017

@author: DanielM
"""

from ouropy.genneuron import GenNeuron
import ouropy.parameters as params


class MossyCell(GenNeuron):
    """Implements the HippCell class with logic inherited from
    ouropy.GenNeuron, a generic neuron class"""
    name = "MossyCell"

    def __init__(self, name=None):
        self.name = name
        self.all_secs = []
        self.dendrites = []

        self.mk_soma(name='soma', diam=20, L=20)
        # Set up the sections with topology
        """NOTE! There is a discrepancy here between the diameter of the third
        section of the dendrite. Publication: 25; Code: 2.5 The number in
        the code seems to be more physiological"""
        self.mk_dendrite(4, dend_name='dend_1',
                         sec_names=['proxd', 'mid1d', 'mid2d', 'dd'],
                         diam=[5.78, 4, 2.5, 1], L=[50, 50, 50, 50],
                         soma_loc=1)
        self.mk_dendrite(4, dend_name='dend_2',
                         sec_names=['proxd', 'mid1d', 'mid2d', 'dd'],
                         diam=[5.78, 4, 2.5, 1], L=[50, 50, 50, 50],
                         soma_loc=1)
        self.mk_dendrite(4, dend_name='dend_3',
                         sec_names=['proxd', 'mid1d', 'mid2d', 'dd'],
                         diam=[5.78, 4, 2.5, 1], L=[50, 50, 50, 50],
                         soma_loc=0)
        self.mk_dendrite(4, dend_name='dend_4',
                         sec_names=['proxd', 'mid1d', 'mid2d', 'dd'],
                         diam=[5.78, 4, 2.5, 1], L=[50, 50, 50, 50],
                         soma_loc=0)

        # Difference between mossycell and mossycell_cat.
        # mossycell_cat imports mossycellparams_cat.txt
        parameters = params.read_parameters('mossycellparams_cat.txt')

        self.insert_mechs(parameters)
