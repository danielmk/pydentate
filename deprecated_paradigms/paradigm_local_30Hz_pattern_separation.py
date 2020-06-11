# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""


from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
import net_tunedrev
from burst_generator_inhomogeneous_poisson import inhom_poiss
import os
import argparse
import scipy.stats as stats

# Handle command line inputs
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
parser.add_argument('-scale',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=500,
                    dest='input_scale')
parser.add_argument('-seed',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=10000,
                    dest='seed')

args = parser.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
input_scale = args.input_scale
seed = args.seed

# Where to search for nrnmech.dll file. Must be adjusted for your machine.
dll_files = [("C:\\Users\\DanielM\\Repos\\models_dentate\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll"),
             "C:\\Users\\daniel\\Repos\\nrnmech.dll",
             ("C:\\Users\\Holger\\danielm\\models_dentate\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll"),
             ("C:\\Users\\Daniel\\repos\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll")]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

np.random.seed(seed)
# Generate a gaussian probability density function for PP->GC mapping.
gauss_gc = stats.norm(loc=1000, scale=input_scale)
gauss_bc = stats.norm(loc=12, scale=(input_scale/2000.0)*24)
pdf_gc = gauss_gc.pdf(np.arange(2000))
pdf_gc = pdf_gc/pdf_gc.sum()
pdf_bc = gauss_bc.pdf(np.arange(24))
pdf_bc = pdf_bc/pdf_bc.sum()
# We hold the pdf constant. To randomize the centroid we reslice the GC indices
GC_indices = np.arange(2000)
start_idc = np.random.randint(0, 1999, size=400)

PP_to_GCs = []
for x in start_idc:
    curr_idc = np.concatenate((GC_indices[x:2000], GC_indices[0:x]))
    PP_to_GCs.append(np.random.choice(curr_idc, size=100, replace=False,
                                      p=pdf_gc))

PP_to_GCs = np.array(PP_to_GCs)
# Generate the PP -> BC mapping as above
BC_indices = np.arange(24)
start_idc = np.array(((start_idc/2000.0)*24), dtype=int)

PP_to_BCs = []
for x in start_idc:
    curr_idc = np.concatenate((BC_indices[x:24], BC_indices[0:x]))
    PP_to_BCs.append(np.random.choice(curr_idc, size=1, replace=False,
                                      p=pdf_bc))

PP_to_BCs = np.array(PP_to_BCs)

# Generate temporal patterns for the 100 PP inputs
np.random.seed(seed)
temporal_patterns = inhom_poiss(rate=30)

dt = 0.1
# Start the runs of the model
for run in runs:
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

    tuned_save_file_name = (str(nw) + "_data_paradigm_local-pattern" +
                            "-separation-30Hz_run_scale_seed_" +
                            str(run).zfill(3) + "_" +
                            str(input_scale).zfill(3) + "_" + str(seed))
    nw.shelve_network(savedir, tuned_save_file_name)

    fig = nw.plot_aps(time=600)
    tuned_fig_file_name = (str(nw) + "_spike-plot_paradigm_local-pattern" +
                           "-separation-30Hz_run_scale_seed_" +
                           str(run).zfill(3) + "_" +
                           str(input_scale).zfill(3) + "_" + str(seed))
    nw.save_ap_fig(fig, savedir, tuned_fig_file_name)
