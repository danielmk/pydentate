# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h
import ouropy
import numpy as np
from granulecell import GranuleCell
from mossycell_cat import MossyCell
from basketcell import BasketCell
from hippcell import HippCell

h.nrn_load_dll("C:\\Users\\Holger\\danielm\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")

class MFStimNetwork(ouropy.gennetwork.GenNetwork):
    
    def __init__(self, seed = None, n_cells = 200, stim_int = 20):
        #Setup cells
        self.mk_population(GranuleCell, 500)
        self.mk_population(MossyCell, 15) 
        self.mk_population(BasketCell, 6)
        self.mk_population(HippCell, 6)
        
        #Set seed for reproducibility
        self.set_numpy_seed(seed)
        
        # Setup recordings
        self.populations[0].record_aps()
        self.populations[1].record_aps()
        self.populations[2].record_aps()
        self.populations[3].record_aps()

        # Set up artificial mf stimulation
        self.iclamped_cells = self.populations[0].current_clamp_rnd(n_cells, amp=1, dur=6, delay=np.arange(100,100+10*stim_int,stim_int))


        self.vclamp_vec = h.Vector()
        vclamp_cell = np.random.choice(self.populations[0].cells, 1, replace = False)[0]
        while vclamp_cell in self.iclamped_cells:
            vclamp_cell = np.random.choice(self.populations[0].cells, 1, replace = False)[0]
        
        vclamp_cell._vclamp()
        self.vclamp_vec.record(vclamp_cell.vclamp._ref_i)
        self.volt_record = vclamp_cell._voltage_recording()

        """
        Call signature of mk_Exp2SynConnection:
        (self, pre_pop, post_pop, target_pool,
         target_segs, divergence, tau1, tau2, e, thr, delay, weight)
        """
        
        # GC -> MC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[1], 3,
                                   'proxd', 1, 0.5,6.2, 0, 10, 1.5, 0.2*10**(-3))
        
        # GC -> BC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[2], 3,
                                   'proxd', 1, 0.3, 0.6, 0, 10, 0.8, 4.7*10**(-3))
        
        # GC -> HC
        self.mk_Exp2SynConnection(self.populations[0], self.populations[3], 5,
                                   'proxd', 3, 0.3, 0.6, 0, 10, 1.5, 0.5*10**(-3))
        
        # MC -> GC
        """self.mk_Exp2SynConnection(self.populations[1], self.populations[0], 350,
                                   'proxd', 200, 1.5, 5.5, 0, 10, 3, 0.3*10**(-3))"""
        
        # MC -> MC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[1], 6,
                                   'proxd', 1, 0.45, 2.2, 0, 10, 2, 0.5*10**(-3))
        
        # MC -> BC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[2], 3,
                                   'proxd', 1, 0.1, 0.1, 1, 10, 3, 0.3*10**(-3))
        
        # MC -> HC
        self.mk_Exp2SynConnection(self.populations[1], self.populations[3], 5,
                                   'midd', 2, 0.9, 3.6, 0, 10, 3,0.2*10**(-3))
        
        # BC -> GC
        #ORIGINAL
        """self.mk_Exp2SynConnection(self.populations[2], self.populations[0], 140,
                                   'soma', 100, 0.26, 5.5, -70, -10, 0.85, 1.6*10**(-3))"""
        
        self.mk_Exp2SynConnection(self.populations[2], self.populations[0], 140,
                                   'soma', 100, 1, 18, -70, -10, 0.85, 1.6*10**(-3))

        # BC -> MC
        self.mk_Exp2SynConnection(self.populations[2], self.populations[1], 7,
                                   'proxd', 3, 0.3, 3.3, -70, -10, 1.5, 1.5*10**(-3))
        
        # BC -> BC
        self.mk_Exp2SynConnection(self.populations[2], self.populations[2], 3,
                                   'proxd', 2, 0.16, 1.8, -70, -10, 0.8, 7.6*10**(-3))
        
        # HC -> GC
        #ORIGINAL
        """self.mk_Exp2SynConnection(self.populations[3], self.populations[0], 260,
                                   'dd', 160, 0.5, 6, -70, 10, 1.6, 0.5*10**(-3))"""
        self.mk_Exp2SynConnection(self.populations[3], self.populations[0], 260,
                                   'dd', 160, 1, 18, -70, 10, 1.6, 0.5*10**(-3))
        
        # HC -> MC
        self.mk_Exp2SynConnection(self.populations[3], self.populations[1], 5,
                                   ['mid1d', 'mid2d'], 4, 0.5, 6, -70, 10, 1, 1.5*10**(-3))
        
        # HC -> BCa
        self.mk_Exp2SynConnection(self.populations[3], self.populations[2], 5,
                                   'ddend', 4, 0.4, 5.8, -70, 10, 1.6, 0.5*10**(-3))

for trial in range(20):
    print("Initializing trial " + str(trial))
    nw = MFStimNetwork(seed=10000+trial, n_cells=200, stim_int=33)

    """Initialization for -2000 to -100"""
    print("Running trial " + str(trial))
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
    h.secondorder = 2
    h.t = 0
    h.dt = 0.01

    """Setup run control for -100 to 1500"""
    h.frecord_init()  # Necessary after changing t to restart the vectors
    while h.t < 500:
        h.fadvance()

    print("Saving trial " + str(trial))
    nw.populations[0].write_aps("GranuleCells_" + str(trial))
    nw.populations[1].write_aps("MossyCells_" + str(trial))
    nw.populations[2].write_aps("BasketCells_" + str(trial))
    nw.populations[3].write_aps("HIPPCells_" + str(trial))

    f_name = "volt_clamp_" + str(trial)

    np.savez(f_name, nw.vclamp_vec.as_numpy())
