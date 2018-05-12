# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h
import numpy as np
import net_tuned
import matplotlib.pyplot as plt
import os

# Locate and load the nrnmech.dll file. Must to be adjusted for your machine.
dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
            "C:\\Users\\daniel\\repos\\nrnmech.dll"]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

# Setup specs for stimulation
n_cells = 60  # Number of cells that are stimulated
stim_pool = 150  # Size of the pool from which stimulated cells are chosen
stim_location = int(2000 / 2.0 - stim_pool / 2.0)
stim_amp = 1
stim_dur = 5
stim_delay = 50
stim_ints = [20]

# Setup specs for measurements
cells_to_measure = np.arange(0, 2000, 50)

save_dir = "C:\\Users\\DanielM\\Repos\\pyDentate\\paradigm_frequency_inhibition_saves_2018-05-09"

for run in range(10):
    for interval in stim_ints:
        # Create a standard networks and add the stimulation
        nw = net_tuned.TunedNetwork(seed=10000+run)
        np.random.seed(10000 + run)
    
        # Make sure we are not stimulating a cell we measure
        stim_cells = np.random.choice(2000,n_cells)
        gcs_to_measure = np.random.choice(2000,20)
        while np.intersect1d(stim_cells, gcs_to_measure).any():
            gcs_to_measure = np.random.choice(2000,20)
    
        for x in range(10):
            nw.populations[0].current_clamp_range(stim_cells,
                                                  amp=stim_amp,
                                                  dur=stim_dur,
                                                  delay=100+interval*x)

        nw.populations[0].SEClamp(gcs_to_measure, dur1 = 100+200+interval*10, rs=1)
        #nw.populations[1].SEClamp([30], dur1 = 100+200+interval*10, rs=1)
        #nw.populations[2].SEClamp([12], dur1 = 100+200+interval*10, rs=1)
        #nw.populations[3].SEClamp([12], dur1 = 100+200+interval*10, rs=1)
        nw.populations[0].voltage_recording(range(0,2000,100))
        nw.populations[1].voltage_recording(range(0,60,2))
        nw.populations[2].voltage_recording(range(24))
        nw.populations[3].voltage_recording(range(24))

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
        while h.t < 100+200+interval*10:
            h.fadvance()
    
        spike_plot = nw.plot_aps(time=100+200+interval*10)
        spike_plot_file_name = "run_" + str(run) + "_spike_plot_intervals_" + str(interval)
        nw.save_ap_fig(spike_plot, directory = save_dir, file_name = spike_plot_file_name + '_' + str(nw))
        data_file_name = "run_" + str(run) + "_data_intervals_"+ str(interval)
        nw.shelve_network(directory = save_dir, file_name = data_file_name + '_' + str(nw))
