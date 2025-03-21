"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
from pydentate import net_tunedrev, neuron_tools
from pydentate.inputs import inhom_poiss
import os
import argparse
import scipy.stats as stats
import platform

# Handle command line inputs
pr = argparse.ArgumentParser(description='Local pattern separation paradigm')
pr.add_argument('-runs',
                nargs=3,
                type=int,
                help='start stop range for the range of runs',
                default=[0, 1, 1],
                dest='runs')
pr.add_argument('-savedir',
                type=str,
                help='complete directory where data is saved',
                default=os.getcwd(),
                dest='savedir')
pr.add_argument('-scale',
                type=int,
                help='standard deviation of gaussian distribution',
                default=1000,
                dest='input_scale')
pr.add_argument('-input_seed',
                type=int,
                help='input_seed',
                default=[10000],
                dest='input_seed')
pr.add_argument('-network_seed',
                type=int,
                help='standard deviation of gaussian distribution',
                default=[10000],
                dest='nw_seed')
pr.add_argument('-input_frequency',
                type=int,
                help='standard deviation of gaussian distribution',
                default=[10],
                dest='input_frequency')

args = pr.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
input_scale = args.input_scale
nw_seed = args.nw_seed
input_seed = args.input_seed
input_frequency = args.input_frequency

# Where to search for nrnmech.dll file. Must be adjusted for your machine.
"""
dirname = os.path.dirname(__file__)
if platform.system() == 'Windows':
    dll_dir = os.path.join(dirname, 'win64', 'nrnmech.dll')
else:
    dll_dir = os.path.join(dirname, 'x86_64', 'libnrnmech.so')    
print("DLL loaded from: " + dll_dir)
# h.nrn_load_dll(dll_dir)
""" 
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Start the runs of the model
for run in runs:
    # Seed the numpy random number generator for replication
    np.random.seed(input_seed[0]+run)

    # Randomly choose target cells for the PP lines
    gauss_gc = stats.norm(loc=1000, scale=input_scale)
    gauss_bc = stats.norm(loc=12, scale=(input_scale/2000.0)*24)
    pdf_gc = gauss_gc.pdf(np.arange(2000))
    pdf_gc = pdf_gc/pdf_gc.sum()
    pdf_bc = gauss_bc.pdf(np.arange(24))
    pdf_bc = pdf_bc/pdf_bc.sum()
    GC_indices = np.arange(2000)
    start_idc = np.random.randint(0, 1999, size=400)

    PP_to_GCs = []
    for x in start_idc:
        curr_idc = np.concatenate((GC_indices[x:2000], GC_indices[0:x]))
        PP_to_GCs.append(np.random.choice(curr_idc, size=100, replace=False,
                                          p=pdf_gc))

    PP_to_GCs = np.array(PP_to_GCs)
    PP_to_GCs = PP_to_GCs[0:24]

    BC_indices = np.arange(24)
    start_idc = np.array(((start_idc/2000.0)*24), dtype=int)

    PP_to_BCs = []
    for x in start_idc:
        curr_idc = np.concatenate((BC_indices[x:24], BC_indices[0:x]))
        PP_to_BCs.append(np.random.choice(curr_idc, size=1, replace=False,
                                          p=pdf_bc))

    PP_to_BCs = np.array(PP_to_BCs)
    PP_to_BCs = PP_to_BCs[0:24]

    # Generate temporal patterns for the 100 PP inputs
    temporal_patterns = inhom_poiss(modulation_rate=input_frequency)
    temporal_patterns[0:24]
    nw = net_tunedrev.TunedNetwork(nw_seed[0], temporal_patterns,
                                   PP_to_GCs,
                                   PP_to_BCs)

    # Attach voltage recordings to all cells
    nw.populations[0].voltage_recording(range(2000))
    nw.populations[1].voltage_recording(range(60))
    nw.populations[2].voltage_recording(range(24))
    nw.populations[3].voltage_recording(range(24))
    # Run the model
    """Initialization for -2000 to -100"""
    print("Running model")
    neuron_tools.run_neuron_simulator()

    tuned_save_file_name = (str(nw) + "-data-paradigm-local-pattern" +
                            "-separation_nw-seed_input-seed_input-frequency_scale_run_" +
                            str(nw_seed[0]) + '_' +
                            str(input_seed[0]) + '_' + 
                            str(input_frequency[0]) + '_' + 
                            str(input_scale).zfill(3) + '_' +
                            str(run).zfill(3) + '_')

    nw.shelve_network(savedir, tuned_save_file_name)

    fig = nw.plot_aps(time=600)
    tuned_fig_file_name = (str(nw) + "_spike-plot_paradigm_local-pattern" +
                           "-separation_run_scale_seed_input-seed_nw-seed_" +
                           str(run).zfill(3) + '_' +
                           str(input_scale).zfill(3) + '_' + str(10000) +
                           str(input_seed) + str(nw_seed))
    nw.save_ap_fig(fig, savedir, tuned_fig_file_name)

