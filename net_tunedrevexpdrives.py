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
gennet = ouropy.gennetwork 

class TunedNetwork(ouropy.gennetwork.GenNetwork):
    """ This model implements the ring model from Santhakumar et al. 2005.
    with some changes as in Yim et al. 2015.
    It features inhibition but omits the MC->GC connection.
    """
    name = "TunedNetworkExpDrives"

    def __init__(self, seed=None, n_gcs=2000, n_mcs=60, n_bcs=24, n_hcs=24,
                 W_pp_gc=1e-3, W_pp_bc=1e-3, n_pp_gc = 20, n_pp_bc = 20, W_gc_bc=2.5e-2, W_gc_hc=2.5e-2,
                 W_bc_gc=1.2e-3, W_hc_gc=6e-3, temporal_patterns=np.array([]), rec_cond=True):
        self.seed=seed
        # Set seed for reproducibility
        if seed:
            self.set_numpy_seed(seed)

        self.init_params = locals()
        self.init_params['self'] = str(self.init_params['self'])
        # Setup cells
        self.mk_population(GranuleCell, n_gcs)
        self.mk_population(MossyCell, n_mcs)
        self.mk_population(BasketCell, n_bcs)
        self.mk_population(HippCell, n_hcs)

        # Setup recordings
        self.populations[0].record_aps()
        self.populations[1].record_aps()
        self.populations[2].record_aps()
        self.populations[3].record_aps()

        t_patterns = np.array(temporal_patterns)

        # PP -> GC
        ouropy.gennetwork.ImplicitConvergentTmgsynConnectionExpProb(self.populations[0], t_patterns, 'midd', n_pp_gc,
                     10, 0, 1, 0, 0, W_pp_gc, rec_cond=rec_cond)

        # PP -> BC
        ouropy.gennetwork.ImplicitConvergentTmgsynConnectionExpProb(self.populations[2], t_patterns, 'ddend', n_pp_bc,
                     6.3, 0, 1, 0, 0, W_pp_bc, rec_cond=rec_cond)

        # GC -> MC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[0], self.populations[1],
                                  6/60, 'proxd',
                                  1, 7.6, 500, 0.1, 0, 0, 10, 1.5, 2e-2, rec_cond=rec_cond)

        # GC -> BC
        #Weight x4, target_pool = 2
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[0], self.populations[2],
                                           4/24, 'proxd',
                                           1, 8.7, 500, 0.1, 0, 0, 10, 0.8, W_gc_bc, rec_cond=rec_cond)

        # GC -> HC
        # Divergence x4; Weight doubled; Connected randomly.
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[0], self.populations[3],
                                           24/24, 'proxd',
                                           1, 8.7, 500, 0.1, 0, 0, 10, 1.5, W_gc_hc, rec_cond=rec_cond)

        # MC -> MC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[1], self. populations[1],
                                           24/60, 'proxd',
                                           3, 2.2, 0, 1, 0, 0, 10, 2, 0.5e-3, rec_cond=rec_cond)

        # MC -> BC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[1], self.populations[2],
                                           6/24, 'proxd',
                                           1, 2, 0, 1, 0, 0, 10, 3, 0.3e-3, rec_cond=rec_cond)

        # MC -> HC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[1], self.populations[3],
                                           10/24, 'midd',
                                           2, 6.2, 0, 1, 0, 0, 10, 3, 0.2e-3, rec_cond=rec_cond)

        # BC -> GC
        # Nr. synapses x3; Weight *1/4; changed from 5.5 to 20 (Hefft & Jonas, 2005)
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[2], self.populations[0],
                                           280/2000, 'soma',
                                           400, 20, 0, 1, 0, -70, 10, 0.85, W_bc_gc, rec_cond=rec_cond)

        # BC -> MC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[2], self.populations[1],
                                           14/60, 'proxd',
                                           3, 3.3, 0, 1, 0, -70, 10, 1.5, 1.5e-3, rec_cond=rec_cond)

        # BC -> BC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[2], self.populations[2],
                                           6/24,'proxd',
                                           2, 1.8, 0,1,0,-70, 10, 0.8, 7.6e-3, rec_cond=rec_cond)

        # HC -> GC
        # Weight x10; Nr synapses x4; changed from 6 to 20 (Hefft & Jonas, 2005)
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[3], self.populations[0],
                                           2000/2000, 'dd',
                                           640, 20, 0, 1, 0, -70, 10, 3.8, W_hc_gc, rec_cond=rec_cond)

        # HC -> MC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[3], self.populations[1],
                                           30/60, ['mid1d', 'mid2d'],
                                           4, 6, 0, 1, 0, -70, 10, 1, 1.5e-3, rec_cond=rec_cond)

        # HC -> BC
        ouropy.gennetwork.DivergentTmgsynConnectionExpProb(self.populations[3], self.populations[2],
                                           12/24, 'ddend',
                                           4, 5.8, 0, 1, 0, -70, 10, 1.6, 0.5e-3, rec_cond=rec_cond)
