# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h
from mfstimnetwork import MFStimNetwork
import numpy as np
import time

for trial in range(20):
    print("Initializing trial " + str(trial))
    nw = MFStimNetwork(seed=10000+trial, n_cells=200, stim_int=20)

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
    while h.t < 400:
        h.fadvance()

    print("Saving trial " + str(trial))
    nw.populations[0].write_aps("GranuleCells_" + str(trial))
    nw.populations[1].write_aps("MossyCells_" + str(trial))
    nw.populations[2].write_aps("BasketCells_" + str(trial))
    nw.populations[3].write_aps("HIPPCells_" + str(trial))

    f_name = "volt_clamp_" + str(trial)

    np.savez(f_name, nw.vclamp_vec.as_numpy())

    stop_save = time.clock()
