# -*- coding: utf-8 -*-
"""
The output files have the following structure
nw-name_run_seed_input-freq_input-scale_bc-decay_hc-decay_bc-delay_hc-delay_gc-weight_bc-weight_hc-weight

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
import net_tunedrevdecaysdelays
from burst_generator_inhomogeneous_poisson import inhom_poiss
import os
import argparse
import scipy.stats as stats

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
pr.add_argument('-seeds',
                nargs=3,
                type=int,
                help='seeds to go through',
                default=[10000, 10001, 1],
                dest='seeds')
pr.add_argument('-input_freq',
                type=int,
                help='frequency of inputs',
                default=10,
                dest='input_freq')
pr.add_argument('-bc_decay',
                type=int,
                help='decay time constant of BC to GC synapse in ms',
                default=20,
                dest='bc_decay')
pr.add_argument('-hc_decay',
                type=int,
                help='decay time constant of HC to GC synapse in ms',
                default=20,
                dest='hc_decay')
pr.add_argument('-bc_delay',
                type=float,
                help='decay time constant of BC to GC synapse in ms',
                default=0.85,
                dest='bc_delay')
pr.add_argument('-hc_delay',
                type=float,
                help='decay time constant of HC to GC synapse in ms',
                default=3.8,
                dest='hc_delay')
pr.add_argument('-gc_weight',
                type=float,
                help='weight of GC to BC and HC synapses',
                default=2.5*10**(-2),
                dest='gc_weight')
pr.add_argument('-bc_weight',
                type=float,
                help='weight of BC to GC synapses',
                default=1.2*10**(-3),
                dest='bc_weight')
pr.add_argument('-hc_weight',
                type=float,
                help='weight of HC to GC synapses',
                default=0.6*10**(-2),
                dest='hc_weight')

args = pr.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
input_scale = args.input_scale
seeds = range(args.seeds[0], args.seeds[1], args.seeds[2])
input_freq = args.input_freq
bc_decay = args.bc_decay
hc_decay = args.hc_decay
bc_delay = args.bc_delay
hc_delay = args.hc_delay
gc_weight = args.gc_weight
bc_weight = args.bc_weight
hc_weight = args.hc_weight

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
print("DLL loaded from: " + dll_dir)
h.nrn_load_dll(dll_dir)

for seed in seeds:
    print("Seed: " + str(seed))
    # Seed the numpy random number generator for replication
    np.random.seed(seed)
    
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
    temporal_patterns = inhom_poiss(rate=10)

    # Start the runs of the model
    for run in runs:
        print("Run: " + str(run))
        nw = net_tunedrevdecaysdelays.TunedNetwork(seed+run,
                                                   temporal_patterns[0+run:24+run],
                                                   PP_to_GCs[0+run:24+run],
                                                   PP_to_BCs[0+run:24+run],
                                                   bc_decay,
                                                   hc_decay,
                                                   bc_delay,
                                                   hc_delay,
                                                   gc_weight,
                                                   bc_weight,
                                                   hc_weight)
        break
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
        h.finitialize(-60)
        h.t = -2000
        h.secondorder = 0
        h.dt = 10
        while h.t < -100:
            h.fadvance()

        h.secondorder = 2
        h.t = 0
        h.dt = 0.1

        """Setup run control for -100 to 600"""
        h.frecord_init()  # Necessary after changing t to restart the vectors
        while h.t < 600:
            h.fadvance()
        print("Done Running")

        save_file_name = (str(nw) + '_' +
                          str(run).zfill(3) + '_' +
                          str(seed).zfill(5) + '_' +
                          str(input_freq) + '_' +
                          str(input_scale).zfill(3) + '_' +
                          str(bc_decay).zfill(3) + '_' + 
                          str(hc_decay).zfill(3) + '_' +
                          str(bc_delay).zfill(6) + '_' + 
                          str(hc_delay).zfill(6) + '_' +
                          str(gc_weight).zfill(6) + '_' + 
                          str(bc_weight).zfill(6) + '_' + 
                          str(hc_weight).zfill(6))
    
        ap_time_stamps = [np.array(x[0]) for x in nw.populations[0].ap_counters]
        ap_binary_array = np.zeros((2000,int(600/dt)), dtype=np.uint8)
        for idx, x in enumerate(ap_time_stamps):
            if x.any():
                spike_idc = np.array(x / dt, dtype=np.int)
                ap_binary_array[idx,spike_idc] = 1
    
        np.savez(savedir + '\\' + 'output_' + save_file_name, ap_binary_array)
        
        gc_inputs = np.zeros((2000,int(600/dt)),dtype=np.uint8)
    
        for idx_pp, pp in enumerate(PP_to_GCs):
            for gc in pp:
                for times in temporal_patterns[idx_pp]:
                    gc_inputs[gc][int(times/dt)] = gc_inputs[gc][int(times/dt)] + 1
        
        np.savez(savedir + '\\' + 'input_' + save_file_name, gc_inputs)
        
        fig = nw.plot_aps(time=600)
        nw.save_ap_fig(fig, savedir, 'figure_' + save_file_name)
        
        if run==0:
            nw.shelve_network(savedir, 'nw-specs_' + save_file_name)