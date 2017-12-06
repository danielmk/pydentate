# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h, gui
from granulecell import GranuleCell
from mossycell_cat import MossyCell
from basketcell import BasketCell
from hippcell import HippCell
import matplotlib.pyplot as plt
import ouropy
import numpy as np

h.nrn_load_dll("C:\Users\DanielM\Repos\models_dentate\dentate_gyrus_Santhakumar2005_and_Yim_patterns\dentategyrusnet2005\\nrnmech.dll")

class PatternNetwork(ouropy.gennetwork.GenNetwork):
    """ This model implements the ring model from Santhakumar et al. 2005.
    It features inhibition but omits the MC->GC connection.
    Perforant Path input is delivered to the first 100 granule cells.
    """

    def __init__(self, temporal_pattern, seed=None, sprouting=0):
        # Setup cells
        # The model is scaled up x4
        """self.mk_population(GranuleCell, 2000)
        self.mk_population(MossyCell, 60)
        self.mk_population(BasketCell, 24)
        self.mk_population(HippCell, 24)"""
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

        self.test_recordings = self.populations[0].voltage_recording(10)

        # Target poools have to be scaled up as well
        # PP -> GC
        self.mk_PerforantPathPoissonStimulation(self.populations[0], temporal_pattern,
                                         200, 'dd',
                                         1.5, 5.5, 0, 2*10**(-2))

        # PP -> MC
        self.mk_PerforantPathPoissonStimulation(self.populations[1], temporal_pattern,
                                         2, 'dd',
                                         1.5, 5.5, 0, 0.5*10**(-2))

        # PP -> BC
        self.mk_PerforantPathPoissonStimulation(self.populations[2], temporal_pattern,
                                         2, 'ddend',
                                         2, 6.3, 0, 1*10**(-2))
        
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
                                  3, 0.45, 2.2, 0, 10, 2, 0.5*10**(-3))

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

        # Sprouting
        self.mk_Exp2SynConnection(self.populations[0], self.populations[0],
                                  25, 'proxd', sprouting,
                                  1.5, 5.5, 0, 10, 0.8, 2*10**(-3))

if __name__ == '__main__':

    np.random.seed(10000)
    # Generate temporal patterns for the 100 PP inputs
    temporal_patterns = np.random.poisson(10, (100, 3)).cumsum(axis = 1)

    # Generate the PP -> GC mapping so that each GC receives inputs from 20/100
    # randomly chosen PP inputs
    innervation_pattern = np.array([np.random.choice(100,20, replace = False) for x in range(500)])
    innervation_pattern = innervation_pattern.swapaxes(0,1)

    PP_to_GCs = []
    for x in range(0,100):
        PP_to_GCs.append(np.argwhere(innervation_pattern == x)[:,1])

    PP_to_GCs = np.array(PP_to_GCs)
    all_targets = np.array([y for x in PP_to_GCs for y in x])

    time_vec = h.Vector()

    nw = PatternNetwork(temporal_patterns[0], seed=10000, sprouting=0)
    time_vec.record(h._ref_t)
    
    """stimvec = h.Vector(temporal_patterns[0])
    nw.connections[0].vecstim.play(stimvec)
    nw.connections[1].vecstim.play(stimvec)
    nw.connections[2].vecstim.play(stimvec)"""

    """Initialization for -2000 to -100"""
    h.cvode.active(0)
    dt = 0.1
    h.steps_per_ms = 1.0/dt
    h.tstop = 1500
    h.finitialize(-60)
    h.t = -2000
    h.secondorder = 0
    h.dt = 10
    while h.t < -100:
        h.fadvance()
        print(h.t)

    h.secondorder = 2
    h.t = 0
    h.dt = 0.1

    """Setup run control for -100 to 1500"""
    h.frecord_init() # Necessary after changing t to restart the vectors

    while h.t < 300:
        h.fadvance()
