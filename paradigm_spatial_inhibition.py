# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h
import ouropy
import numpy as np
import net_tuned
import time

# Office PC
h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")
#Home PC
#h.nrn_load_dll("C:\\Users\\daniel\\repos\\nrnmech.dll")

# Setup specs for stimulation
n_cells = 100 # Number of cells that are stimulated
stim_pool = 150 # Size of the pool from which stimulated cells are chosen
stim_location = int(2000 / 2.0 - stim_pool / 2.0)
stim_amp = 1
stim_dur = 6
stim_delay = 100

# Setup specs for measurements
cells_to_measure = range(0,2000, 5)

for run in range(1):
    # Create a standard networks and add the stimulation
    nw = net_tuned.TunedNetwork(seed = 10000 + run)
    stim_cells = np.random.choice(range(stim_location,stim_location+stim_pool), n_cells, replace = False)
    nw.populations[0].current_clamp_range(stim_cells,
                                               amp = stim_amp,
                                               dur = stim_dur,
                                               delay = stim_delay)
    
    nw.populations[0].SEClamp(cells_to_measure)
    nw.populations[0].voltage_recording(cells_to_measure)

    """measured_cells = []
    for x in cells_to_measure:
        ctrl_vclamp_cell = nww.populations[0].cells[x]
        ctrl_vclamp_cell._SEClamp(dur1=200, amp1=0, rs=0.001)()
        ctrl_i_vec = h.Vector()
        ctrl_i_vec.record(ctrl_vclamp_cell.vclamp._ref_i)
        ctrl_volt_rec = ctrl_vclamp_cell._voltage_recording()"""
    
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

    """nw.populations[0].write_aps("GranuleCells_trial_" + str(trial) + "stim_pos_" + str(stim_pos))
    nw.populations[1].write_aps("MossyCells_" + str(trial) + "stim_pos_" + str(stim_pos))
    nw.populations[2].write_aps("BasketCells_" + str(trial) + "stim_pos_" + str(stim_pos))
    nw.populations[3].write_aps("HIPPCells_" + str(trial) + "stim_pos_" + str(stim_pos))
    
    f_name = "volt_clamp_trial_" + str(trial) + "stim_pos_" + str(stim_pos)
    
    np.savez(f_name, nw.vclamp_vec.as_numpy())"""
