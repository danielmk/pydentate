# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h
import numpy as np
import standard_net
import time

h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")

# Setup specs of the paradigm
n_cells = 800  # Number of cells that are stimulated
stim_amp = 1
stim_dur = 6
stim_delay = 100
idx_measured_mc = 12
stim_periods = [1000, 100, 33, 20] # in ms
#stim_periods = [1000]
reps = 1  # Number of repetitions for the networks

for period in stim_periods:
    stim_times = np.arange(100, 100+10*period, period)
    for rep in range(reps):
        # Create a standard networks and add the stimulation
        ctrl_nw = standard_net.StandardNetwork(seed = 10000 + rep)
        #targets = np.random.choice(2000, n_cells, replace = False)
        ctrl_nw.populations[0].mk_current_clamp(n_cells,
                                               amp = stim_amp,
                                               dur = stim_dur,
                                               delays = stim_times)
        ctrl_vclamp_cell = ctrl_nw.populations[1].cells[12]
        ctrl_vclamp_cell._vclamp(dur1=100+10*period+500, amp1=-70, rs=0.001)
        ctrl_i_vec = h.Vector()
        ctrl_i_vec.record(ctrl_vclamp_cell.vclamp._ref_i)
        ctrl_volt_rec = ctrl_vclamp_cell._voltage_recording()

        """Initialization for -2000 to -100"""
        #print("Running rep " + str(rep))
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
        while h.t < stim_times[9]+500:
            h.fadvance()

        print("Saving rep " + str(rep))
        ctrl_nw.populations[0].write_aps("GranuleCells_period_" + str(period) + "_rep_" + str(rep))
        ctrl_nw.populations[1].write_aps("MossyCells_period_" + str(period) + "_rep_" + str(rep))
        ctrl_nw.populations[2].write_aps("BasketCells_period_" + str(period) + "_rep_" + str(rep))
        ctrl_nw.populations[3].write_aps("HIPPCells_period_" + str(period) + "_rep_" + str(rep))
        
        f_name = "volt_clamp_period_" + str(period) + "_rep_" + str(rep)
        
        np.savez(f_name, ctrl_i_vec.as_numpy())
        
        stop_save = time.clock()
