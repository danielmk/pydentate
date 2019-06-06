# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace 
import numpy as np
import net_tunedrevexpdrives
from input_generator import inhom_poiss
import os
import argparse

# Handle command line inputs with argparse
parser = argparse.ArgumentParser(description='Pattern separation paradigm')
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

args = parser.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir

# Where to search for nrnmech.dll file. Must be adjusted for your machine.
dll_files = [("/home/daniel/repos/pyDentate/mechs_7-6_linux/x86_64/.libs/libnrnmech.so")]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

# Generate temporal patterns for the 100 PP inputs
np.random.seed(10000)
temporal_patterns = inhom_poiss()

# Start the runs of the model
for run in runs:
    for idx in range(run+24,temporal_patterns.shape[0]):
        temporal_patterns[idx] = np.array([])

    nw = net_tunedrevexpdrives.TunedNetwork(seed=10000, temporal_patterns=temporal_patterns)

    # Run the model
    """Initialization for -2000 to -100"""
    h.cvode.active(0)
    dt = 0.1
    h.steps_per_ms = 1.0/dt
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

    while h.t < 600:
        h.fadvance()
    print("Done Running")

    tuned_save_file_name = str(nw) + '_run_' + str(run)
    nw.shelve_network(savedir, tuned_save_file_name)

    fig = nw.plot_aps(time=600)
    tuned_fig_file_name = str(nw) + '_spike_plot_run_' + str(run)
    nw.save_ap_fig(fig, savedir, tuned_fig_file_name)
