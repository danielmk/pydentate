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
import time

h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")

# Setup specs of the paradigm
n_cells = 160 # Number of cells that are stimulated
stim_pool = 200 # Size of the pool from which stimulated cells are chosen
stim_amp = 1
stim_dur = 6
stim_delay = 100
idx_measured_cell = 500
offset = 50 # Step size between stim positions
stim_start = range(0, 2001 - stim_pool, offset)
reps = 1 # Number of repetitions for the networks

for trial in stim_start:
    for rep in range(reps):
        # Create a standard networks and add the stimulation
        ctrl_nw = standard_net.StandardNetwork(seed = 10000 + rep)
        target_pool = range(trial, trial + 200)
        targets = np.random.choice(target_pool, n_cells, replace = False)
        targets = np.delete(targets,np.argwhere(targets == idx_measured_cell))
        ctrl_nw.populations[0].current_clamp_range(targets,
                                                   amp = stim_amp,
                                                   dur = stim_dur,
                                                   delay = stim_delay)
        ctrl_vclamp_cell = ctrl_nw.populations[0].cells[500]
        ctrl_vclamp_cell._vclamp()
        ctrl_i_vec = h.Vector()
        ctrl_i_vec.record(ctrl_vclamp_cell.vclamp._ref_i)
        ctrl_volt_rec = ctrl_vclamp_cell._voltage_recording()
        
        """Initialization for -2000 to -100"""
        #print("Running trial " + str(trial))
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
        while h.t < 200:
            h.fadvance()

        print("Saving trial " + str(trial))
        nw.populations[0].write_aps("GranuleCells_trial_" + str(trial) + "stim_pos_" + str(stim_pos))
        nw.populations[1].write_aps("MossyCells_" + str(trial) + "stim_pos_" + str(stim_pos))
        nw.populations[2].write_aps("BasketCells_" + str(trial) + "stim_pos_" + str(stim_pos))
        nw.populations[3].write_aps("HIPPCells_" + str(trial) + "stim_pos_" + str(stim_pos))
        
        f_name = "volt_clamp_trial_" + str(trial) + "stim_pos_" + str(stim_pos)
        
        np.savez(f_name, nw.vclamp_vec.as_numpy())
        
        stop_save = time.clock()
