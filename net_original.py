# -*- coding: utf-8 -*-
"""
This module implements the class StandardNetwork.
StandardNetwork creates a ring network as defined in Santhakumar et al. 2005
with some changes as in Yim et al. 2015.
See StandardNetwork docstring for details.
Created on Tue Nov 28 13:01:38 2017

@author: DanielM
"""

from neuron import h, gui
import ouropy
import matplotlib.pyplot as plt
import numpy as np
from granulecell import GranuleCell
from mossycell_cat import MossyCell
from basketcell import BasketCell
from hippcell import HippCell

h.nrn_load_dll("C:\Users\DanielM\Repos\models_dentate\dentate_gyrus_Santhakumar2005_and_Yim_patterns\dentategyrusnet2005\\nrnmech.dll")

class StandardNetworkOriginal(ouropy.gennetwork.GenNetwork):
    """ This model implements the ring model from Santhakumar et al. 2005.
    with some changes as in Yim et al. 2015.
    It features inhibition but omits the MC->GC connection.
    """

    def __init__(self, seed=None, temporal_patterns=np.array([]),
                 spatial_patterns_gcs=np.array([]),
                 spatial_patterns_bcs=np.array([]), sprouting=0):
        # Setup cells
        self.mk_population(GranuleCell, 500)
        self.mk_population(MossyCell, 15)
        self.mk_population(BasketCell, 6)
        self.mk_population(HippCell, 6)

        # Set seed for reproducibility
        if seed:
            self.set_numpy_seed(seed)

        # Setup recordings
        self.populations[0].record_aps()
        self.populations[1].record_aps()
        self.populations[2].record_aps()
        self.populations[3].record_aps()
        
        temporal_patterns = np.atleast_2d(temporal_patterns)
        if spatial_patterns_gcs.any() and temporal_patterns.any():
            spatial_patterns_gcs = np.atleast_2d(spatial_patterns_gcs)
            for pat in range(len(spatial_patterns_gcs)):

                ouropy.gennetwork.PerforantPathPoissonStimulation(self.populations[0],
                                                           temporal_patterns[pat],
                                                           spatial_patterns_gcs[pat],
                                                           'dd',
                                                           1.5, 5.5, 0, 2*10**(-2))

        if spatial_patterns_bcs.any() and temporal_patterns.any():
            spatial_patterns_bcs = np.atleast_2d(spatial_patterns_bcs)
            for pat in range(len(spatial_patterns_bcs)):
                # PP -> BC
                ouropy.gennetwork.PerforantPathPoissonStimulation(self.populations[2],
                                                           temporal_patterns[pat],
                                                           spatial_patterns_bcs[pat],
                                                           'ddend',
                                                           2, 6.3, 0, 1*10**(-2))

        # Sprouting
        ouropy.gennetwork.Exp2SynConnection(self.populations[0], self.populations[0],
                                  100, 'proxd', sprouting,
                                  1.5, 5.5, 0, 10, 0.8, 2*10**(-3))
        
        # GC -> MC
        ouropy.gennetwork.Exp2SynConnection(self.populations[0], self.populations[1],
                                  3, 'proxd',
                                  1, 0.5,6.2, 0, 10, 1.5, 0.2*10**(-3))
    
        # GC -> BC
        ouropy.gennetwork.Exp2SynConnection(self.populations[0], self.populations[2],
                                  3, 'proxd',
                                  1, 0.3, 0.6, 0, 10, 0.8, 4.7*10**(-3))

        # GC -> HC
        ouropy.gennetwork.Exp2SynConnection(self.populations[0], self.populations[3],
                                  5, 'proxd',
                                  3, 0.3, 0.6, 0, 10, 1.5, 0.5*10**(-3))

        # MC -> MC
        ouropy.gennetwork.Exp2SynConnection(self.populations[1], self.populations[1],
                                  6, 'proxd',
                                  3, 0.45, 2.2, 0, 10, 2, 0.5*10**(-3))

        # MC -> BC
        ouropy.gennetwork.Exp2SynConnection(self.populations[1], self.populations[2],
                                  3, 'proxd',
                                  1, 0.1, 0.1, 1, 10, 3, 0.3*10**(-3))

        # MC -> HC
        ouropy.gennetwork.Exp2SynConnection(self.populations[1], self.populations[3],
                                  5, 'midd',
                                  2, 0.9, 3.6, 0, 10, 3,0.2*10**(-3))

        # BC -> GC
        #ORIGINAL
        ouropy.gennetwork.Exp2SynConnection(self.populations[2], self.populations[0],
                                     140, 'soma',
                                     100, 0.26, 5.5, -70, -10, 0.85, 1.6*10**(-3))

        # BC -> MC
        ouropy.gennetwork.Exp2SynConnection(self.populations[2], self.populations[1],
                                  7, 'proxd',
                                  3, 0.3, 3.3, -70, -10, 1.5, 1.5*10**(-3))

        # BC -> BC
        ouropy.gennetwork.Exp2SynConnection(self.populations[2], self.populations[2],
                                  3, 'proxd',
                                  2, 0.16, 1.8, -70, -10, 0.8, 7.6*10**(-3))

        # HC -> GC
        #ORIGINAL
        ouropy.gennetwork.Exp2SynConnection(self.populations[3], self.populations[0],
                                     260, 'dd',
                                     160, 0.5, 6, -70, 10, 1.6, 0.5*10**(-3))

        # HC -> MC
        ouropy.gennetwork.Exp2SynConnection(self.populations[3], self.populations[1],
                                  5, ['mid1d', 'mid2d'],
                                  4, 0.5, 6, -70, 10, 1, 1.5*10**(-3))

        # HC -> BC
        ouropy.gennetwork.Exp2SynConnection(self.populations[3], self.populations[2],
                                  5, 'ddend',
                                  4, 0.4, 5.8, -70, 10, 1.6, 0.5*10**(-3))
