#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 14:48:40 2019

@author: barisckuru
"""


import numpy as np
from mossycell import MossyCell
from basketcell import BasketCell
import os
from neuron import h, gui

# load modified dll files
dll = "/home/can/repos/pyDentate/mechs_7-6_linux/x86_64/.libs/libnrnmech.so"
dll_files = [(dll)]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

'''SETUP/SIMULATE THE MOSSY CELL AND SYNAPTIC PROCESS

simulate() function in “tmgexp2_simulator” simulates stimulation of a
specific type of cell in NEURON and connects a specific type of synapse onto
that cell. Then it stimulates the cell according to given tau facilitation, tau
recovery and frequency; optionally u0 and sampling rate. Recording are set up.
After simulation it finds the peaks of the conductance (G) response. Finally,
it returns G response and peaks for each stimulation in simulation '''


def simulate(tau_facil, tau_rec, freq, u0=None, sampling=None):
    # u0 & sampling rate are optional varibles
    # MossyCell created
    mycell = MossyCell()
    # Attach tmgexp2syn to MossyCell, 
    # secs[0] is the "proximal dendrite" as stated in Santhakumar 2005 paper
    sec = mycell.dendrites[0].secs[0]
    syn = h.tmgexp2syn(sec(0.5))
    # set the tau rise & decay, and dendrite for specific synapse
    # tau_1(tau rise) & tau_2(tau decay)
    # tau_1=0.5, tau_2=6.2 values gc to mc from the Santhakumar 2005 paper
    # tau_1=0.3, tau_2=0.6 values gc to in from the Santhakumar 2005 paper
    syn.tau_1 = 0.3
    syn.tau_2 = 0.6
    syn.tau_facil = tau_facil
    syn.tau_rec = tau_rec
    stim_period_ms = 1000/freq
    # 0.078 is optimized u0 value for 2 different data sets
    if u0 is None:
        u0 = 0.078
    syn.u0 = u0
    stim = h.NetStim()
    stim.interval = stim_period_ms
    # set the sampling rate here, 0.05 is the sampling rate of recorded data
    if sampling is None:
        sampling = 0.5
    sample_period_ms = sampling
    # set starting time of stimuli in the simulation here
    # stim_start to be used further
    stim_start = 1000
    stim.start = stim_start
    # set the number of stimuli
    stim.number = 10
    ncstim = h.NetCon(stim, syn)
    # ncstim.delay can be used to set the delay
    # set the stimuli weight
    ncstim.weight[0] = 0.001
    # variable simulation duration calculated for each freq
    simdur_init = 16000
    simdur = simdur_init/freq
    # Setup recording from MossyCell soma with h.SEClamp()
    c = h.SEClamp(mycell.soma(0.5))
    c.dur1 = simdur+stim_start
    c.amp1 = -70.42
    c.rs = 0.1
    v_vec = h.Vector()
    t_vec = h.Vector()
    i_vec = h.Vector()
    g_syn_vec = h.Vector()
    t_vec.record(h._ref_t)
    v_vec.record(mycell.soma(0.5)._ref_v)
    i_vec.record(c._ref_i)
    g_syn_vec.record(syn._ref_g)
    # Run Simulation
    h.cvode.active(0)
    h.dt = sample_period_ms
    h.steps_per_ms = 1.0/h.dt
    h.finitialize(-70.42)
    h.secondorder = 0
    while h.t < simdur+stim_start:
        h.fadvance()
    # Return output arrays
    ivec = np.array(i_vec)
    vvec = np.array(v_vec)
    tvec = np.array(t_vec)
    # Cut the beginning, we needed it for active dyn to settle down
    cut = 500  # ms
    start = int(cut/sample_period_ms)
    end_cut = (len(tvec)-int(100/sample_period_ms))
    ivec = ivec[start:end_cut]
    vvec = vvec[start:end_cut]
    tvec = tvec[start:end_cut]
    stim_start = stim_start-cut
    gvec = -np.divide(ivec, vvec)
    # Offset
    gvec = gvec-(np.average(gvec[0:cut]))
    # Normalization of simulation data by the peak of the first signal
    stim_start_dtp = int(stim_start/sample_period_ms)
    stim_stop_dtp = int(stim_start_dtp+stim_period_ms/sample_period_ms)
    sim_first_max = max(gvec[stim_start_dtp:stim_stop_dtp])
    norm_sim = gvec/sim_first_max

    'PEAKS'
    # period in data points
    period = stim_period_ms/sample_period_ms
    # first & last stimuli in data points
    first_stim = stim_start/sample_period_ms
    last_stim = (stim_start+int(stim_period_ms)*(stim.number+1))/sample_period_ms
    # peaks indices aranged
    peaks_idc = np.arange(first_stim, last_stim, period, dtype=int)
    # split the each stimuli in simulation
    split_sim = np.split(norm_sim, peaks_idc)
    split_sim = split_sim[1:-1:]
    split_sim = np.array(split_sim)
    # find the peaks
    peaks_sim = np.amax(split_sim, axis=1)
    return peaks_sim, norm_sim
