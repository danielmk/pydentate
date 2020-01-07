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
import matplotlib.pyplot as plt
import numpy as np
from granulecell import GranuleCell
from mossycell_cat import MossyCell
from basketcell import BasketCell
from hippcell import HippCell
from ouropy.gennetwork import GenNetwork, ConnDivergent, ImplicitConvergentTmgsynConnectionExpProb

class TunedNetwork(GenNetwork):
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
        ImplicitConvergentTmgsynConnectionExpProb(self.populations[0], t_patterns, 'midd', n_pp_gc,
                     1.5, 10, 2.22590184e+13, 0.1, 3.12395254e+02, 8.15226996e-01, 0, W_pp_gc, rec_cond=rec_cond)

        # PP -> BC
        ImplicitConvergentTmgsynConnectionExpProb(self.populations[2], t_patterns, 'ddend', n_pp_bc,
                     2, 6.3, 0, 1, 0, 0, 0, W_pp_bc, rec_cond=rec_cond)

        # GC -> MC
        gc_mc_syn = {'tau_1': 0.5,
                     'tau_2': 8.7,
                     'tau_facil': 3.73987038e+03,
                     'U': 0.1,
                     'tau_rec': 3.96206731e+01,
                     'u0': 5.94309891e-02,
                     'e': 0,}
        gc_mc_net = {'threshold': 10,
                    'delay': 1.5,
                    'weight': 2e-2}

        ConnDivergent(self.populations[0], self.populations[1],
                                  6/60, 'proxd', 1, h.tmgexp2syn, gc_mc_syn,  gc_mc_net, rec_cond=rec_cond)

        # GC -> BC
        #Weight x4, target_pool = 2
        gc_bc_syn = {'tau_1': 0.3,
                     'tau_2': 8.7,
                     'tau_facil': 3.52671493e+03,
                     'U': 0.1,
                     'tau_rec': 1.33249553e+01,
                     'u0': 5.77352002e-02,
                     'e': 0,}
        gc_bc_net = {'threshold': 10,
                    'delay': 0.8,
                    'weight':  W_gc_bc}
        ConnDivergent(self.populations[0], self.populations[2],
                                           4/24, 'proxd', 1, h.tmgexp2syn, gc_bc_syn, gc_bc_net, rec_cond=rec_cond)

        # GC -> HC
        # Divergence x4; Weight doubled; Connected randomly.
        gc_hc_syn = {'tau_1': 0.3,
                     'tau_2': 8.7,
                     'tau_facil': 3.52671493e+03,
                     'U': 0.1,
                     'tau_rec': 1.33249553e+01,
                     'u0': 5.77352002e-02,
                     'e': 0,}
        gc_hc_net = {'threshold': 10,
                    'delay': 1.5,
                    'weight':  W_gc_hc}
        ConnDivergent(self.populations[0], self.populations[3],
                                           24/24, 'proxd',
                                           1, h.tmgexp2syn, gc_hc_syn, gc_hc_net, rec_cond=rec_cond)

        # MC -> GC
        mc_gc_syn = {'tau_1': 1.5,
                     'tau_2': 5.5,
                     'tau_facil': 0,
                     'U': 0.1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': 0,}
        mc_gc_net = {'threshold': 10,
                    'delay': 3,
                    'weight':  3e-4}
        ConnDivergent(self.populations[1], self. populations[0],
                                           42/60, 'proxd',200, h.tmgexp2syn, mc_gc_syn, mc_gc_net, rec_cond=rec_cond)

        # MC -> MC
        mc_mc_syn = {'tau_1': 0.45,
                     'tau_2': 2.2,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': 0,}
        mc_mc_net = {'threshold': 10,
                    'delay': 2,
                    'weight':  5e-4}
        ConnDivergent(self.populations[1], self. populations[1],
                                           24/60, 'proxd',
                                           3, h.tmgexp2syn, mc_mc_syn, mc_mc_net, rec_cond=rec_cond)

        # MC -> BC
        mc_bc_syn = {'tau_1': 0.9,
                     'tau_2': 2,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': 0,}
        mc_bc_net = {'threshold': 10,
                    'delay': 3,
                    'weight':  3e-4}
        ConnDivergent(self.populations[1], self.populations[2],
                                           6/24, 'proxd',
                                           1, h.tmgexp2syn, mc_bc_syn, mc_bc_net, rec_cond=rec_cond)

        # MC -> HC
        mc_hc_syn = {'tau_1': 0.9,
                     'tau_2': 6.2,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': 0,}
        mc_hc_net = {'threshold': 10,
                    'delay': 3,
                    'weight': 2e-4}
        ConnDivergent(self.populations[1], self.populations[3],
                                           10/24, 'midd',
                                           2, h.tmgexp2syn, mc_hc_syn, mc_hc_net, rec_cond=rec_cond)

        # BC -> GC
        # Nr. synapses x3; Weight *1/4; changed from 5.5 to 20 (Hefft & Jonas, 2005)
        bc_gc_syn = {'tau_1': 0.26,
                     'tau_2': 20,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': -70,}
        bc_gc_net = {'threshold': 10,
                    'delay': 0.85,
                    'weight': W_bc_gc}
        ConnDivergent(self.populations[2], self.populations[0],
                                           280/2000, 'soma',
                                           400, h.tmgexp2syn, bc_gc_syn, bc_gc_net, rec_cond=rec_cond)

        # BC -> MC
        bc_mc_syn = {'tau_1': 0.3,
                     'tau_2': 3.3,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': -70,}
        bc_mc_net = {'threshold': 10,
                    'delay': 1.5,
                    'weight':1.5e-3}
        ConnDivergent(self.populations[2], self.populations[1],
                                           14/60, 'proxd',
                                           3, h.tmgexp2syn, bc_mc_syn, bc_mc_net, rec_cond=rec_cond)

        # BC -> BC
        bc_bc_syn = {'tau_1': 0.16,
                     'tau_2': 1.8,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': -70,}
        bc_bc_net = {'threshold': 10,
                    'delay': 0.8,
                    'weight':1.5e-3}
        ConnDivergent(self.populations[2], self.populations[2],
                                           6/24,'proxd',
                                           2, h.tmgexp2syn, bc_bc_syn, bc_bc_net, rec_cond=rec_cond)

        # HC -> GC
        # Weight x10; Nr synapses x4; changed from 6 to 20 (Hefft & Jonas, 2005)
        hc_gc_syn = {'tau_1': 0.5,
                     'tau_2': 20,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': -70,}
        hc_gc_net = {'threshold': 10,
                    'delay': 3.8,
                    'weight': W_hc_gc}
        ConnDivergent(self.populations[3], self.populations[0],
                                           2000/2000, 'dd',
                                           640, h.tmgexp2syn, hc_gc_syn, hc_gc_net, rec_cond=rec_cond)

        # HC -> MC
        hc_mc_syn = {'tau_1': 0.5,
                     'tau_2': 6,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': -70,}
        hc_mc_net = {'threshold': 10,
                    'delay': 1,
                    'weight': 1.5e-3}
        ConnDivergent(self.populations[3], self.populations[1],
                                           30/60, ['mid1d', 'mid2d'],
                                           4, h.tmgexp2syn, hc_mc_syn, gc_mc_net, rec_cond=rec_cond)

        # HC -> BC
        hc_bc_syn = {'tau_1': 0.4,
                     'tau_2': 5.8,
                     'tau_facil': 0,
                     'U': 1,
                     'tau_rec': 0,
                     'u0': 0,
                     'e': -70,}
        hc_bc_net = {'threshold': 10,
                    'delay': 1.6,
                    'weight': 5e-4}
        ConnDivergent(self.populations[3], self.populations[2],
                                           12/24, 'ddend',
                                           4, h.tmgexp2syn, hc_bc_syn, hc_bc_net, rec_cond=rec_cond)
