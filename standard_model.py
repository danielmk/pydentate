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
from ouropy.gennetwork import GenNetwork
import numpy as np

h.nrn_load_dll("C:\Users\DanielM\Repos\models_dentate\dentate_gyrus_Santhakumar2005\dentategyrusnet2005\\nrnmech.dll")

#Setup the cells of the network
celltypes = [GranuleCell, MossyCell, BasketCell, HippCell]
cellnums = [500,15,6,6]
my_nw = GenNetwork(celltypes, cellnums)

# Setup recordings
gc_ap_counter = my_nw.populations[0].record_aps()
mc_ap_counter = my_nw.populations[1].record_aps()
bc_ap_counter = my_nw.populations[2].record_aps()
hc_ap_counter = my_nw.populations[3].record_aps()
time = h.Vector()
time.record(h._ref_t)

gc_vectors = []
for x in range(0,100,20):
    gc_vectors.append(my_nw.populations[0].cells[x]._voltage_recording())

# Set up current clamp to simulate frequency stim
# 33 Hz

# Use IClamp to simulate mossy fiber stimulation
iclamped_cells = my_nw.populations[0].current_clamp_rnd(200, amp=1, dur=8, delay=np.arange(100,300,30))

vclamp_vec = h.Vector()
vclamp_cell = np.random.choice(my_nw.populations[0].cells, 1, replace = False)[0]
while vclamp_cell in iclamped_cells:
    vclamp_cell = np.random.choice(my_nw.populations[0].cells, 1, replace = False)[0]

vclamp_cell._vclamp()
vclamp_vec.record(vclamp_cell.vclamp._ref_i)
volt_record = vclamp_cell._voltage_recording()

"""
Call signature of PerforantPathStimulation:
(self, post_pop, n_targets, target_segs,
                 tau1, tau2, e, thr, delay, weight)
"""

# PP -> GC
"""my_nw.mk_PerforantPathStimulation(my_nw.populations[0], np.arange(100), 'dd',
                 1.5, 5.5, 0, 10, 3, 2*10**(-2))"""

# PP -> BC
"""my_nw.mk_PerforantPathStimulation(my_nw.populations[2], 2, 'ddend',
                 2, 6.3, 0, 10, 3, 1*10**(-2))"""
                 
# PP -> MC
"""my_nw.mk_PerforantPathStimulation(my_nw.populations[1], 2, 'dd',
                 1.5, 5.5, 0, 10, 3, 0.5*10**(-2))"""

# Sprouting
"""my_nw.mk_Sprouting(my_nw.populations[0],10,'proxd',
                           1.5, 5.5, 0, 10, 0.8, 2*10**(-2))"""

"""
Call signature of mk_Exp2SynConnection:
(self, pre_pop, post_pop, target_pool,
 target_segs, divergence, tau1, tau2, e, thr, delay, weight)
"""

# GC -> MC
my_nw.mk_Exp2SynConnection(my_nw.populations[0], my_nw.populations[1], 4,
                           'proxd', 1, 0.5,6.2, 0, 10, 1.5, 0.2*10**(-3))

# GC -> BC
my_nw.mk_Exp2SynConnection(my_nw.populations[0], my_nw.populations[2], 4,
                           'proxd', 1, 0.3, 0.6, 0, 10, 0.8, 4.7*10**(-3))

# GC -> HC
my_nw.mk_Exp2SynConnection(my_nw.populations[0], my_nw.populations[3], 5,
                           'proxd', 3, 0.3, 0.6, 0, 10, 1.5, 0.5*10**(-3))

# MC -> GC
"""my_nw.mk_Exp2SynConnection(my_nw.populations[1], my_nw.populations[0], 4,
                           'proxd', 1, 1.5, 5.5, 0, 10, 3, 0.3*10**(-3))"""

# MC -> MC
my_nw.mk_Exp2SynConnection(my_nw.populations[1], my_nw.populations[1], 6,
                           'proxd', 1, 0.45, 2.2, 0, 10, 2, 0.5*10**(-3))

# MC -> BC
my_nw.mk_Exp2SynConnection(my_nw.populations[1], my_nw.populations[2], 3,
                           'proxd', 1, 0.1, 0.1, 1, 10, 3, 0.3*10**(-3))

# MC -> HC
my_nw.mk_Exp2SynConnection(my_nw.populations[1], my_nw.populations[3], 5,
                           'midd', 2, 0.9, 3.6, 0, 10, 3,0.2*10**(-3))

# BC -> GC
my_nw.mk_Exp2SynConnection(my_nw.populations[2], my_nw.populations[0], 140,
                           'soma', 100, 0.26, 5.5, -70, -10, 0.85, 1.6*10**(-3))

# BC -> MC
my_nw.mk_Exp2SynConnection(my_nw.populations[2], my_nw.populations[1], 7,
                           'proxd', 3, 0.3, 3.3, -70, -10, 1.5, 1.5*10**(-3))

# BC -> BC
my_nw.mk_Exp2SynConnection(my_nw.populations[2], my_nw.populations[1], 3,
                           'proxd', 2, 0.16, 1.8, -70, -10, 0.8, 7.6*10**(-3))

# HC -> GC
my_nw.mk_Exp2SynConnection(my_nw.populations[3], my_nw.populations[0], 260,
                           'dd', 160, 0.5, 6, -70, 10, 1.6, 0.5*10**(-3))

# HC -> MC
my_nw.mk_Exp2SynConnection(my_nw.populations[3], my_nw.populations[1], 5,
                           ['mid1d', 'mid2d'], 4, 0.5, 6, -70, 10, 1, 1.5*10**(-3))

# HC -> BC
my_nw.mk_Exp2SynConnection(my_nw.populations[3], my_nw.populations[2], 5,
                           'ddend', 4, 0.4, 5.8, -70, 10, 1.6, 0.5*10**(-3))

#Measure a synapse
syn_vec = h.Vector()
syn_vec.record(my_nw.connections[9].synapses[3][0]._ref_i)

"""
for x in range(0,500,100):
    my_nw.populations[0].cells[x]._current_clamp_soma(amp=0.3, dur=20, delay=300)
"""

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
h.dt = 0.01

"""Setup run control for -100 to 1500"""
h.frecord_init() # Necessary after changing t to restart the vectors
while h.t < 600:
    h.fadvance()
