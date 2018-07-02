# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h
import numpy as np
import net_tunedrev
from burst_generator_inhomogeneous_poisson import inhom_poiss, hom_poiss
import os
import argparse
import scipy.stats as stats

# Handle command line inputs with argparse
parser = argparse.ArgumentParser(description='Local pattern separation paradigm')
parser.add_argument('-runs',
                    nargs=3,
                    type=int,
                    help='start stop range for the range of runs',
                    default=[0, 1, 1],
                    dest='runs')
parser.add_argument('-savedir',
                    type=str,
                    help='complete directory where data is saved',
                    default=os.getcwd(),
                    dest='savedir')
parser.add_argument('-rate',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=200,
                    dest='input_rate')
parser.add_argument('-seed',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=10000,
                    dest='seed')

args = parser.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
input_rate = args.input_rate
seed = args.seed

# Locate a nrnmech.dll file that has the mechanisms required by the network
# On your own machine you have to add the path to your own file to the list dll_files
dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
            "C:\\Users\\daniel\\Repos\\nrnmech.dll",
            "C:\\Users\\Holger\\danielm\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
            "C:\\Users\\Daniel\\repos\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll"]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

# Generate temporal patterns for the 100 PP inputs
np.random.seed(seed)
temporal_patterns = hom_poiss(input_rate)

innervation_pattern_gc = np.array([np.random.choice(400,20, replace = False) for x in range(2000)])
innervation_pattern_gc = innervation_pattern_gc.swapaxes(0,1)

PP_to_GCs = []
for x in range(0,400):
    PP_to_GCs.append(np.argwhere(innervation_pattern_gc == x)[:,1])

PP_to_GCs = np.array(PP_to_GCs)

# Generate the PP -> BC mapping as above
innervation_pattern_bc = np.array([np.random.choice(400,20, replace = False) for x in range(24)])
innervation_pattern_bc = innervation_pattern_bc.swapaxes(0,1)

PP_to_BCs = []
for x in range(0,400):
    PP_to_BCs.append(np.argwhere(innervation_pattern_bc == x)[:,1])

PP_to_BCs = np.array(PP_to_BCs)
all_targets = np.array([y for x in PP_to_GCs for y in x])

# Start the runs of the model
for run in runs:
    print("Starting Setup for run " + str(run))
    nw = net_tunedrev.TunedNetwork(seed, temporal_patterns[0+run:24+run],
                                PP_to_GCs[0+run:24+run],
                                PP_to_BCs[0+run:24+run])

    # Attach voltage recordings to all cells
    nw.populations[0].voltage_recording(range(2000))
    nw.populations[1].voltage_recording(range(60))
    nw.populations[2].voltage_recording(range(24))
    nw.populations[3].voltage_recording(range(24))

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

    h.secondorder = 2
    h.t = 0
    h.dt = 0.1
    
    print("Starting run " + str(run))

    """Setup run control for -100 to 1500"""
    h.frecord_init() # Necessary after changing t to restart the vectors

    while h.t < 600:
        h.fadvance()

    print("Done Running " + str(run))

    tuned_save_file_name = str(nw) + '_data_paradigm_rate-pattern-separation_run_scale_seed_' + str(run).zfill(3) + '_' + str(input_rate).zfill(3) + '_' + str(seed)
    nw.shelve_network(savedir, tuned_save_file_name)

    fig = nw.plot_aps(time=600)
    tuned_fig_file_name = str(nw) + '_spike-plot_paradigm_rate-pattern-separation_run_scale_seed_' + str(run).zfill(3) + '_' + str(input_rate).zfill(3) + '_' + str(seed)
    nw.save_ap_fig(fig, savedir, tuned_fig_file_name)
