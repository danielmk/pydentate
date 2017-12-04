# -*- coding: utf-8 -*-
"""
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

h.nrn_load_dll("C:\Users\Holger\danielm\models_dentate\dentate_gyrus_Santhakumar2005_and_Yim_patterns\dentategyrusnet2005\\nrnmech.dll")

class LightStimNetwork(ouropy.gennetwork.GenNetwork):
    """ This model implements the ring model from Santhakumar et al. 2005.
    It features inhibition but omits the MC->GC connection.
    Perforant Path input is delivered to the first 100 granule cells.
    """
    def __init__(self, seed=None, n_cells=40, target_pool = range(40)):
        # Setup cells
        self.mk_population(GranuleCell, 2000)
        self.mk_population(MossyCell, 60)
        self.mk_population(BasketCell, 24)
        self.mk_population(HippCell, 24)

        # Set seed for reproducibility
        if seed:
            self.set_numpy_seed(seed)

        # Setup recordings
        self.populations[0].record_aps()
        self.populations[1].record_aps()
        self.populations[2].record_aps()
        self.populations[3].record_aps()
        
        self.vclamp_vec = h.Vector()
        vclamp_cell = self.populations[0].cells[500]
        vclamp_cell._vclamp()
        self.vclamp_vec.record(vclamp_cell.vclamp._ref_i)
        self.volt_record = vclamp_cell._voltage_recording()
        
        target_pop = []
        for idx in target_pool:
            target_pop.append(self.populations[0].cells[idx])
        
        clamp_pop = np.random.choice(target_pop, n_cells, replace = False)
        while vclamp_cell in clamp_pop:
            clamp_pop = np.random.choice(target_pop, n_cells, replace = False)
        
        for cell in clamp_pop:
            cell._current_clamp_soma(amp=1, dur=6, delay=100)
        
        
        """
        Call signature of mk_Exp2SynConnection:
        (self, pre_pop, post_pop, target_pool,
         target_segs, divergence, tau1, tau2, e, thr, delay, weight)
        """
        
        # GC -> MC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[1],
                                  12, 'proxd',
                                  1, 0.5,6.2, 0, 10, 1.5, 0.2*10**(-3))
    
        # GC -> BC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[2],
                                  12, 'proxd',
                                  1, 0.3, 0.6, 0, 10, 0.8, 4.7*10**(-3))

        # GC -> HC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[3],
                                  20, 'proxd',
                                  3, 0.3, 0.6, 0, 10, 1.5, 0.5*10**(-3))

        # MC -> GC
        """self.mk_Exp2SynConnection(self.populations[1], self.populations[0],
                                     350, 'proxd',
                                     200, 1.5, 5.5, 0, 10, 3, 0.3*10**(-3))"""

        # MC -> MC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[1],
                                  24, 'proxd',
                                  3, 0.45, 2.2, 0, 10, 2, 0.5*10**(-3))

        # MC -> BC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[2],
                                  12, 'proxd',
                                  1, 0.1, 0.1, 1, 10, 3, 0.3*10**(-3))

        # MC -> HC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[3],
                                  20, 'midd',
                                  2, 0.9, 3.6, 0, 10, 3,0.2*10**(-3))
        
        # BC -> GC
        #ORIGINAL
        """self.mk_Exp2SynConnection(self.populations[2], self.populations[0],
                                     140, 'soma',
                                     100, 0.26, 5.5, -70, -10, 0.85, 1.6*10**(-3))"""
        
        self.mk_Exp2SynConnection(self.populations[2], self.populations[0],
                                  560, 'soma',
                                  100, 0.26, 5.5, -70, -10, 0.85, 1.6*10**(-3))
        
        
        # BC -> MC
        self.mk_Exp2SynConnection(self.populations[2], self.populations[1],
                                  28, 'proxd',
                                  3, 0.3, 3.3, -70, -10, 1.5, 1.5*10**(-3))
        
        # BC -> BC
        self.mk_Exp2SynConnection(self.populations[2], self.populations[2],
                                  12, 'proxd',
                                  2, 0.16, 1.8, -70, -10, 0.8, 7.6*10**(-3))
        
        # HC -> GC
        #ORIGINAL
        """self.mk_Exp2SynConnection(self.populations[3], self.populations[0],
                                     260, 'dd',
                                     160, 0.5, 6, -70, 10, 1.6, 0.5*10**(-3))"""
        self.mk_Exp2SynConnection(self.populations[3], self.populations[0],
                                  1040, 'dd',
                                  160, 0.5, 6, -70, 10, 1.6, 0.5*10**(-3))
        
        # HC -> MC
        self.mk_Exp2SynConnection(self.populations[3], self.populations[1],
                                  20, ['mid1d', 'mid2d'],
                                  4, 0.5, 6, -70, 10, 1, 1.5*10**(-3))
        
        # HC -> BCa
        self.mk_Exp2SynConnection(self.populations[3], self.populations[2],
                                  20, 'ddend',
                                  4, 0.4, 5.8, -70, 10, 1.6, 0.5*10**(-3))

if __name__ == '__main__':

    stimulations = 20
    for trial in range(stimulations):
        for rep in range(3):
            target_pool = range(250*trial - 150 * trial, 250*(trial+1) - 150 * (trial+1))

            nw = LightStimNetwork(seed=10000+rep, n_cells=60, target_pool=target_pool)

            h.cvode.active(0)
            dt = 0.01
            h.steps_per_ms = 1.0/dt
            h.tstop = 1500
            h.finitialize(-60)
            h.t = -2000
            h.secondorder = 0
            h.dt = 10
            while h.t < -100:
                h.fadvance()
                #print(h.t)
            h.secondorder = 2
            h.t = 0
            h.dt = 0.01

            """Setup run control for -100 to 1500"""
            h.frecord_init() # Necessary after changing t to restart the vectors
            while h.t < 200:
                h.fadvance()

            print("Saving trial " + str(trial))
            nw.populations[0].write_aps("GranuleCells_trial_" + str(trial) + "_rep_" + str(rep))
            nw.populations[1].write_aps("MossyCells_" + str(trial) + "_rep_" + str(rep))
            nw.populations[2].write_aps("BasketCells_" + str(trial) + "_rep_" + str(rep))
            nw.populations[3].write_aps("HIPPCells_" + str(trial) + "_rep_" + str(rep))

            f_name = "volt_clamp_trial_" + str(trial) + "_rep_" +str(rep)

            np.savez(f_name, nw.vclamp_vec.as_numpy())