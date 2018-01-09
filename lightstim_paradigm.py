# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h
import ouropy
import numpy as np
import standard_net

h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")

# Setup specs of the paradigm
n_cells = 160 # Number of cells that are stimulated
stim_pool = 200 # Size of the pool from which stimulated cells are chosen
stim_amp = 1
stim_dur = 6
stim_delay = 100
idx_measured_cell = 500

for x in range 
# Create a standard networks
nw = standard_net.StandardNetwork()

"""Initialization for -2000 to -100"""
#print("Running trial " + str(trial))
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
h.secondorder = 2
h.t = 0
h.dt = 0.1

"""Setup run control for -100 to 1500"""
h.frecord_init()  # Necessary after changing t to restart the vectors
while h.t < 200:
    h.fadvance()

"""print("Saving trial " + str(trial))
nw.populations[0].write_aps("GranuleCells_trial_" + str(trial) + "stim_pos_" + str(stim_pos))
nw.populations[1].write_aps("MossyCells_" + str(trial) + "stim_pos_" + str(stim_pos))
nw.populations[2].write_aps("BasketCells_" + str(trial) + "stim_pos_" + str(stim_pos))
nw.populations[3].write_aps("HIPPCells_" + str(trial) + "stim_pos_" + str(stim_pos))

f_name = "volt_clamp_trial_" + str(trial) + "stim_pos_" + str(stim_pos)

np.savez(f_name, nw.vclamp_vec.as_numpy())

stop_save = time.clock()
"""