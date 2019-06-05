# -*- coding: utf-8 -*-
"""
The output files have the following structure
nw-name_run_seed_input-freq_input-scale_bc-decay_hc-decay_bc-delay_hc-delay_gc-weight_bc-weight_hc-weight

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
import net_tunedrevenhancedff
from burst_generator_inhomogeneous_poisson import inhom_poiss
import os
import argparse
import scipy.stats as stats
import time
from analysis_main import time_stamps_to_signal
tsts = time_stamps_to_signal

# Handle command line inputs
pr = argparse.ArgumentParser(description='Local pattern separation paradigm')
pr.add_argument('-runs',
                nargs=3,
                type=int,
                help='start stop range for the range of runs',
                default=[0,1,1],
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
                help='network seeds to go through',
                default=[10000, 10001, 1],
                dest='seeds')
pr.add_argument('-input_freq',
                type=int,
                help='frequency of inputs',
                default=10,
                dest='input_freq')
pr.add_argument('-bc_decay',
                type=float,
                help='decay time constant of BC to GC synapse in ms',
                default=1,
                dest='bc_decay')
pr.add_argument('-hc_decay',
                type=float,
                help='decay time constant of HC to GC synapse in ms',
                default=1,
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
runs = range(args.runs[0],args.runs[1],args.runs[2])
save_dir = args.savedir
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
dll_files = [("/home/daniel/repos/pyDentate/mechs/x86_64/.libs/libnrnmech.so"),
             ("C:\\Users\\Daniel\\repos\\pyDentate\\mechs_7-5\\nrnmech.dll"),
             ("C:\\Users\\DanielM\\Repos\\models_dentate\\"
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
    print("Network seed: " + str(seed))
    # Seed the numpy random number generator for replication
    input_seed = seed+1
    np.random.seed(input_seed)

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
    np.random.seed(input_seed)
    temporal_patterns = inhom_poiss(rate=input_freq)
    
    pp_ts = []
    gc_ts = []
    mc_ts = []
    bc_ts = []
    hc_ts = []
    
    pp_to_gc_list = []
    gc_to_hc_list = []
    gc_to_bc_list = []
    hc_to_gc_list = []
    bc_to_gc_list = []

    seed_input_vectors = []
    seed_output_vectors = []
    seed_input_time_stamps = []
    seed_output_time_stamps = []

    # Start the runs of the model
    for run in runs:
        start_proc_t = time.clock()
        print("Run: " + str(run) + ". Total time: " + str(start_proc_t))
        nw = net_tunedrevenhancedff.TunedNetwork(seed,
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
        end_proc_t = time.clock()
        print("Done Running at " + str(end_proc_t) + " after " + str((end_proc_t - start_proc_t)/60) + " minutes")

        output_pop_vector = np.array([x[1].n for x in nw.populations[0].ap_counters])
        seed_output_vectors.append(output_pop_vector)
        input_pop_vector = np.zeros(temporal_patterns.shape[0])
        input_pop_vector[0+run:24+run] = [x.shape[0] for x in temporal_patterns[0+run:24+run]]
        seed_input_vectors.append(input_pop_vector)


        save_fig_name = (str(nw) + '_' +
                          str(run).zfill(3) + '_' +
                          str(seed).zfill(5) + '_' +
                          str(input_seed).zfill(5) + '_' +
                          str(input_freq) + '_' +
                          str(input_scale).zfill(3) + '_' +
                          str(bc_decay).zfill(3) + '_' +
                          str(hc_decay).zfill(3) + '_' +
                          str(bc_delay).zfill(6) + '_' + 
                          str(hc_delay).zfill(6) + '_' +
                          str(gc_weight).zfill(6) + '_' +
                          str(bc_weight).zfill(6) + '_' +
                          str(hc_weight).zfill(6))

        fig = nw.plot_aps(time=600)
        nw.save_ap_fig(fig, save_dir, 'figure_' + save_fig_name)
        pp_lines = np.empty(400, dtype = np.object)
        pp_lines[0+run:24+run] = temporal_patterns[0+run:24+run]
        # pp_ts.append(pp_lines)

        curr_pp_ts = np.array(tsts(pp_lines, dt_signal=0.1, t_start=0, t_stop=600), dtype = np.bool)
        curr_gc_ts = np.array(tsts(nw.populations[0].get_properties()['ap_time_stamps'], dt_signal=0.1, t_start=0, t_stop=600), dtype = np.bool)
        curr_mc_ts = np.array(tsts(nw.populations[1].get_properties()['ap_time_stamps'], dt_signal=0.1, t_start=0, t_stop=600), dtype = np.bool)
        curr_hc_ts = np.array(tsts(nw.populations[2].get_properties()['ap_time_stamps'], dt_signal=0.1, t_start=0, t_stop=600), dtype = np.bool)
        curr_bc_ts = np.array(tsts(nw.populations[3].get_properties()['ap_time_stamps'], dt_signal=0.1, t_start=0, t_stop=600), dtype = np.bool)

        pp_ts.append(curr_pp_ts)
        gc_ts.append(curr_gc_ts)
        mc_ts.append(curr_mc_ts)
        bc_ts.append(curr_bc_ts)
        hc_ts.append(curr_hc_ts)
        
        pp_lines = nw.populations[0].connections[0:24]
        summed_pp_lines = [np.array(x.conductances).sum(axis=0).sum(axis=0) for x in pp_lines]
        pp_to_gc = np.array(summed_pp_lines).sum(axis=0)
        pp_to_gc_list.append(pp_to_gc)
        
        gc_to_hc_list.append(np.array(nw.populations[0].connections[26].conductances).sum(axis=0).sum(axis=0))
        gc_to_bc_list.append(np.array(nw.populations[0].connections[25].conductances).sum(axis=0).sum(axis=0))
        bc_to_gc_list.append(np.array(nw.populations[0].connections[27].conductances).sum(axis=0).sum(axis=0))
        hc_to_gc_list.append(np.array(nw.populations[0].connections[28].conductances).sum(axis=0).sum(axis=0))
        
        
        
    save_data_name = (str(nw) + '_' +
                      str(seed).zfill(5) + '_' +
                      str(input_seed).zfill(5) + '_' +
                      str(input_freq).zfill(3) + '_' +
                      str(input_scale).zfill(3) + '_' +
                      str(bc_decay).zfill(3) + '_' +
                      str(hc_decay).zfill(3) + '_' +
                      str(bc_delay).zfill(6) + '_' + 
                      str(hc_delay).zfill(6) + '_' +
                      str(gc_weight).zfill(6) + '_' +
                      str(bc_weight).zfill(6) + '_' +
                      str(hc_weight).zfill(6))

    np.savez(save_dir + "time-stamps_" + save_data_name,
             pp_ts = np.array(pp_ts),
             gc_ts = np.array(gc_ts),
             mc_ts = np.array(mc_ts),
             bc_ts = np.array(bc_ts),
             hc_ts = np.array(hc_ts))

    np.savez(save_dir + "conductances_" + save_data_name,
             pp_to_gc = np.array(pp_to_gc_list),
             gc_to_hc = np.array(gc_to_hc_list),
             gc_to_bc = np.array(gc_to_bc_list),
             bc_to_gc = np.array(bc_to_gc_list),
             hc_to_gc = np.array(hc_to_gc_list))

    # Calculate the correlation matrices
    n_runs = len(seed_input_vectors)
    input_corrs = []
    output_corrs = []
    for row_idx in range(n_runs):
        for col_idx in range(1+row_idx,n_runs):
            curr_input_corr = stats.pearsonr(seed_input_vectors[row_idx], seed_input_vectors[col_idx])[0]
            curr_output_corr = stats.pearsonr(seed_output_vectors[row_idx], seed_output_vectors[col_idx])[0]
            input_corrs.append(curr_input_corr)
            output_corrs.append(curr_output_corr)

    save_file_name = (str(nw) + '_' +
                      str(seed).zfill(5) + '_' +
                      str(input_seed).zfill(5) + '_' +
                      str(input_freq).zfill(3) + '_' +
                      str(input_scale).zfill(3) + '_' +
                      str(bc_decay).zfill(3) + '_' +
                      str(hc_decay).zfill(3) + '_' +
                      str(bc_delay).zfill(6) + '_' + 
                      str(hc_delay).zfill(6) + '_' +
                      str(gc_weight).zfill(6) + '_' +
                      str(bc_weight).zfill(6) + '_' +
                      str(hc_weight).zfill(6) + '_')

    np.savez(save_dir + "corrs_" + save_file_name, input_corrs = np.array(input_corrs), output_corrs = np.array(output_corrs))
