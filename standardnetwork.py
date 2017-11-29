# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 12:50:14 2017

@author: DanielM
"""

import ouropy
from neuron import h, gui
from granulecell import GranuleCell
from mossycell_cat import MossyCell
from basketcell import BasketCell
from hippcell import HippCell
import numpy as np

h.nrn_load_dll("C:\Users\DanielM\Repos\models_dentate\dentate_gyrus_Santhakumar2005_and_Yim_patterns\dentategyrusnet2005\\nrnmech.dll")

class StandardNetwork(ouropy.gennetwork.GenNetwork):
    """ This model implements the ring model from Santhakumar et al. 2005.
    It features inhibition but omits the MC->GC connection.
    Perforant Path input is delivered to the first 100 granule cells.
    """
    def __init__(self, seed=None, sprouting=0):
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

        # Setup the PP stimulator
        self.pp_stim = h.NetStim()
        self.pp_stim.number = 1
        self.pp_stim.start = 5

        # PP -> GC
        self.mk_PerforantPathStimulation(self.pp_stim, self.populations[0],
                                         np.arange(100), 'dd',
                                         1.5, 5.5, 0, 10, 3, 2*10**(-2))

        # PP -> MC
        self.mk_PerforantPathStimulation(self.pp_stim, self.populations[1],
                                         2, 'dd',
                                         1.5, 5.5, 0, 10, 3, 0.5*10**(-2))

        # PP -> BC
        self.mk_PerforantPathStimulation(self.pp_stim, self.populations[2],
                                         2, 'ddend',
                                         2, 6.3, 0, 10, 3, 1*10**(-2))

        # Sprouting
        self.mk_Exp2SynConnection(self.populations[0], self.populations[0],
                                  100, 'proxd', sprouting,
                                  1.5, 5.5, 0, 10, 0.8, 2*10**(-3))
        
        """
        Call signature of mk_Exp2SynConnection:
        (self, pre_pop, post_pop, target_pool,
         target_segs, divergence, tau1, tau2, e, thr, delay, weight)
        """
        
        # GC -> MC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[1],
                                  3, 'proxd',
                                  1, 0.5,6.2, 0, 10, 1.5, 0.2*10**(-3))
    
        # GC -> BC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[2],
                                  3, 'proxd',
                                  1, 0.3, 0.6, 0, 10, 0.8, 4.7*10**(-3))

        # GC -> HC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[3],
                                  5, 'proxd',
                                  3, 0.3, 0.6, 0, 10, 1.5, 0.5*10**(-3))

        # MC -> GC
        """self.mk_Exp2SynConnection(self.populations[1], self.populations[0],
                                     350, 'proxd',
                                     200, 1.5, 5.5, 0, 10, 3, 0.3*10**(-3))"""

        # MC -> MC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[1],
                                  6, 'proxd',
                                  1, 0.45, 2.2, 0, 10, 2, 0.5*10**(-3))

        # MC -> BC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[2],
                                  3, 'proxd',
                                  1, 0.1, 0.1, 1, 10, 3, 0.3*10**(-3))

        # MC -> HC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[3],
                                  5, 'midd',
                                  2, 0.9, 3.6, 0, 10, 3,0.2*10**(-3))
        
        # BC -> GC
        #ORIGINAL
        """self.mk_Exp2SynConnection(self.populations[2], self.populations[0],
                                     140, 'soma',
                                     100, 0.26, 5.5, -70, -10, 0.85, 1.6*10**(-3))"""
        
        self.mk_Exp2SynConnection(self.populations[2], self.populations[0],
                                  140, 'soma',
                                  100, 0.26, 20, -70, -10, 0.85, 1.6*10**(-3))
        
        
        # BC -> MC
        self.mk_Exp2SynConnection(self.populations[2], self.populations[1],
                                  7, 'proxd',
                                  3, 0.3, 3.3, -70, -10, 1.5, 1.5*10**(-3))
        
        # BC -> BC
        self.mk_Exp2SynConnection(self.populations[2], self.populations[2],
                                  3, 'proxd',
                                  2, 0.16, 1.8, -70, -10, 0.8, 7.6*10**(-3))
        
        # HC -> GC
        #ORIGINAL
        """self.mk_Exp2SynConnection(self.populations[3], self.populations[0],
                                     260, 'dd',
                                     160, 0.5, 6, -70, 10, 1.6, 0.5*10**(-3))"""
        self.mk_Exp2SynConnection(self.populations[3], self.populations[0],
                                  260, 'dd',
                                  160, 0.5, 20, -70, 10, 1.6, 0.5*10**(-3))
        
        # HC -> MC
        self.mk_Exp2SynConnection(self.populations[3], self.populations[1],
                                  5, ['mid1d', 'mid2d'],
                                  4, 0.5, 6, -70, 10, 1, 1.5*10**(-3))
        
        # HC -> BCa
        self.mk_Exp2SynConnection(self.populations[3], self.populations[2],
                                  5, 'ddend',
                                  4, 0.4, 5.8, -70, 10, 1.6, 0.5*10**(-3))