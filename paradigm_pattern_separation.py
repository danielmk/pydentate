# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h
import ouropy
import numpy as np
import net_tuned
import net_nonfacilitating
import net_global
#import net_tuned_10ECInputs
import os
from burst_generator_inhomogeneous_poisson import inhom_poiss
import sys

print(sys.argv)

dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
            "C:\\Users\\daniel\\Repos\\nrnmech.dll"]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

# Generate temporal patterns for the 100 PP inputs
np.random.seed(10000)
temporal_patterns = inhom_poiss()
# Original from Yim
#np.random.seed(10000)
#temporal_patterns = np.random.poisson(10, (400, 3)).cumsum(axis = 1)

# Generate the PP -> GC mapping so that each GC receives inputs from 20/400
# randomly chosen PP inputs
innervation_pattern_gc = np.array([np.random.choice(400,20, replace = False) for x in range(2000)])
innervation_pattern_gc = innervation_pattern_gc.swapaxes(0,1)

PP_to_GCs = []
for x in range(0,400):
    PP_to_GCs.append(np.argwhere(innervation_pattern_gc == x)[:,1])

PP_to_GCs = np.array(PP_to_GCs)

innervation_pattern_bc = np.array([np.random.choice(400,20, replace = False) for x in range(24)])
innervation_pattern_bc = innervation_pattern_bc.swapaxes(0,1)

PP_to_BCs = []
for x in range(0,400):
    PP_to_BCs.append(np.argwhere(innervation_pattern_bc == x)[:,1])

PP_to_BCs = np.array(PP_to_BCs)
all_targets = np.array([y for x in PP_to_GCs for y in x])

save_dir = "C:\\Users\\DanielM\\Repos\\pyDentate\\paradigm_pattern-separation_saves_2018-05-12_patterns"
    
runs = range(0,20)

if len(sys.argv) > 2:
    runs = range(sys.argv[1],sys.argv[2])
for run in runs:
    nw_tuned = net_tuned.TunedNetwork(10000+run, temporal_patterns[0+run:24+run], PP_to_GCs[0+run:24+run], PP_to_BCs[0+run:24+run], sprouting=0)

    # Attach voltage recordings to all cells
    nw_tuned.populations[0].voltage_recording(range(2000))
    nw_tuned.populations[1].voltage_recording(range(60))
    nw_tuned.populations[2].voltage_recording(range(24))
    nw_tuned.populations[3].voltage_recording(range(24))
    # Run the model
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
    
    while h.t < 1000:
        h.fadvance()
    print("Done Running")

    tuned_save_file_name = str(nw_tuned) + '_run_' + str(run)
    nw_tuned.shelve_network(save_dir, tuned_save_file_name)

    fig = nw_tuned.plot_aps(time=1000)
    tuned_fig_file_name = str(nw_tuned) + '_spike_plot_run_' + str(run)
    nw_tuned.save_ap_fig(fig, save_dir, tuned_fig_file_name)
