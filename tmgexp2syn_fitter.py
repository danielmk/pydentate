# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 21:57:24 2019

@author: Daniel
"""

import numpy as np
from mossycell_cat import MossyCell
from basketcell import BasketCell
from granulecell import GranuleCell
import os
from neuron import h, gui
from scipy.optimize import minimize
import time
import pdb

# load modified dll files
dll_files = [
    "/home/can/repos/pyDentate/mechs_7-6_linux/x86_64/.libs/libnrnmech.so",
    "C:\\Users\\Daniel\\repos\\pyDentate\\mechs_7-6_win\\nrnmech.dll",
]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

def simulate(
    tau_facil,
    tau_rec,
    stim_freq,
    celltype,
    tau_1=0.3,
    tau_2=0.6,
    u0=0.078,
    sampling_interval=0.5,
    v_clamp=-70.42,
    sec='proxd'
):
    """Simulates a frequency stimulation of the tmgexp2syn meachnism onto a
    neuron and measure the conductance at the soma.

    Parameters
    ----------
    tau_facil : numeric
        facilitation time constant of the tmgexp2syn mechanism
        
    tau_rec : numeric
        recovery time constant of the tmgexp2syn mechanism
        
    stim_freq : numeric
        stimulation frequency in Hz
        
    celltype : ouropy.GenNeuron
        A class that inherits from ouropy.GenNeuron
        
    tau_1 : numeric
        Synaptic rise time
        
    tau_2 : numeric
        Synaptic decay time
        
    u0 : float
        The u0 parameters as specified in tmgexp2syn
        
    sampling_interval : numeric
        The sampling interval in ms
        
    v_clamp : numeric
        Holding potential of the
        
    Returns
    -------
    peaks_norm : ndarray
        Peak conductances normalized to first peak
    
    garr_norm : ndarray
        Conductance array normalized to first peak
    """
    # Resolve simulation hyperparameters
    left_pad, right_pad = 50, 100  # in ms
    simdur = left_pad + right_pad + 10 * (1000 / stim_freq)
    # Create cell
    mycell = celltype()
    # Create synapse
    syn = h.tmgexp2syn(mycell.get_segs_by_name(sec)[0](0.5))
    syn.tau_1 = tau_1
    syn.tau_2 = tau_2
    syn.tau_facil = tau_facil
    syn.tau_rec = tau_rec
    syn.u0 = u0
    syn.e = 0
    # Create stimulation
    stim_period_ms = 1000 / stim_freq
    stim = h.NetStim()
    stim.interval = stim_period_ms
    stim.start = left_pad
    stim.number = 10
    nc = h.NetCon(stim, syn)
    nc.weight[0] = 0.001
    # Create voltage clamp
    c = h.SEClamp(mycell.soma(0.5))
    c.dur1 = simdur
    c.amp1 = v_clamp
    c.rs = 0.1
    # Create recording vector
    ivec = h.Vector()
    ivec.record(c._ref_i)
    # Start simulation
    sample_period_ms = sampling_interval
    h.cvode.active(0)
    h.finitialize(-70.42)
    h.t = -2000
    h.secondorder = 0
    h.dt = 10
    while h.t < -100:
        h.fadvance()
    h.secondorder = 0
    h.t = 0
    h.dt = sample_period_ms
    h.steps_per_ms = 1.0 / h.dt
    h.frecord_init()  # Necessary after changing t to restart the vectors
    while h.t < simdur:
        h.fadvance()
    # Get Peaks from output

    iarr = np.array(ivec)
    garr = (iarr - iarr[0:50].mean()) / v_clamp
    stim_dtp = int(np.around(stim.interval/sampling_interval))
    pad_dtp = int(left_pad/sampling_interval)
    t_stim = np.arange(pad_dtp,
                       stim_dtp * 11 + pad_dtp,
                       stim_dtp)
    #dtp_stim = np.array(t_stim / sampling_interval, dtype=np.int)
    dtp_stim = np.array(t_stim, dtype=np.int)
    gvec_split = np.array(np.split(garr, dtp_stim)[1:-1])

    peaks = gvec_split.max(axis=1)
    peaks_norm = peaks / peaks[0]
    garr_norm = garr / peaks[0]

    return peaks_norm, garr_norm


def peakfinder(
    signal, sem, stim, peak_win=(0.002, 0.016), peaks="negative", sampling_rate=20000
):
    """Calculates the peaks normalized to the first.

    Parameters
    ----------
    signal : ndarray
        The signal containing the waveforms.
    sem : ndarray
        Same shape as signal, contains the standard error of mean
    stim : ndarray
        Contains the stimulation points.
    peak_win : sequence of numeric
        Time window after stimulation where peak is expected in s.
    peaks : str, optional
        Direction of the peaks. The default is 'negative'.

    Returns
    -------
    peaks : ndarray
        The normalized peaks.
    peaks_sem_norm : ndarray
        Standard error of means of the peaks normalized to first peak.
    """
    if peaks == "negative":
        signal = -signal
    stim_norm = stim / stim.max()
    stim_thr = (stim_norm.max() - stim_norm.min()) / 2
    stim_idc = np.argwhere(np.diff(stim_norm > stim_thr))[::2, 0]
    peak_win_start = peak_win[0] * sampling_rate
    peak_win_end = peak_win[1] * sampling_rate
    peak_idc = np.concatenate(
        (stim_idc + peak_win_start, stim_idc + peak_win_start + peak_win_end)
    )
    peak_idc = peak_idc.astype(np.int)
    peak_idc.sort()
    signal_split = np.array(np.split(signal, peak_idc)[1:-1:2])
    sem_split = np.array(np.split(sem, peak_idc)[1:-1:2])
    peaks = signal_split.max(axis=1)
    peaks_norm = peaks / peaks[0]
    peak_idc = signal_split.argmax(axis=1)
    peaks_sem = sem_split[range(sem_split.shape[0]), peak_idc]
    peaks_sem_norm = peaks_sem / peaks[0]

    return peaks_norm, peaks_sem_norm


def loss(x, *args):
    """Loss function for use with scipy.optimize.minimize. x is a 1d-array
    that presents the variables that are to be optimized and *args contains all
    other parameters needed to specify the function.
    
    Parameters
    ----------
    x : 1d-array
        x[0], tau_facil. x[1], tau_rec. x[2], u0.
    *args : tuple
        args[0] : vectorized real peaks, array of numerics
        args[1] : celltype class
        args[2] : tau_1 numeric
        args[3] : tau_2 numeric
        args[4] : sampling_interval numeric
        args[5] : v_clamp numeric
        args[6] : sec, h.Section
        args[7] : freqs, sequence of numerics
    """
    # Run the simuations
    sim_peaks = []
    for freq in args[7]:
        sim = simulate(
            x[0],
            x[1],
            freq,
            args[1],
            tau_1=args[2],
            tau_2=args[3],
            u0=x[2],
            sampling_interval=args[4],
            v_clamp=args[5],
            sec=args[6]
        )
        sim_peaks.append(sim[0])


    sim_peaks = np.array(sim_peaks).flatten()

    mean_squared_error = np.square(sim_peaks - args[0]).mean()

    return mean_squared_error


if __name__ == "__main__":
    data_path = (
        "C:\\Users\\Daniel\\Dropbox\\02_MEC Project\\003_Antidromic "
        "electrically evoked EPSCs in hilar cells, voltage clamp\\"
    )
    gc_to_in = np.load(data_path + "gc_to_in_data_full.npz")
    gc_to_mc = np.load(data_path + "gc_to_mc_data_full.npz")
    pp_to_gc = np.load(data_path + "pp_to_gc_data_full.npz")
    """
    # GC to IN FITTING
    gc_to_in_vec = np.concatenate(
        (gc_to_in['peaks_norm_01'], gc_to_in['peaks_norm_10'], gc_to_in['peaks_norm_30'], gc_to_in['peaks_norm_50'])
    )
    gc_to_in_args = (gc_to_in_vec, BasketCell, 0.3, 0.6, 0.5, -70.42, 'proxd', [1,10,30,50])
    gc_to_in_x0 = [3000, 10, 0.1]
    #loss(gc_to_in_x0, gc_to_in_args)
    res_gc_in = minimize(loss, gc_to_in_x0, gc_to_in_args, method='Nelder-Mead')

    # GC to MC FITTING
    gc_to_mc_vec = np.concatenate(
        (gc_to_mc['peaks_norm_01'], gc_to_mc['peaks_norm_10'], gc_to_mc['peaks_norm_30'], gc_to_mc['peaks_norm_50'])
    )
    
    # tau_1=0.5, tau_2=6.2 values gc to mc from the Santhakumar 2005 paper
    gc_to_mc_args = (gc_to_mc_vec, MossyCell, 0.5, 6.2, 0.5, -70.42, 'proxd', [1,10,30,50])
    gc_to_mc_x0 = [3000, 25, 0.1]
    #loss(gc_to_in_x0, gc_to_in_args)
    res_gc_mc = minimize(loss, gc_to_mc_x0, gc_to_mc_args, method='Nelder-Mead')
    """
    pp_to_gc_vec = np.concatenate(
        (pp_to_gc['peaks_norm_05'], pp_to_gc['peaks_norm_10'])
    )

    pp_to_gc_args = (pp_to_gc_vec, GranuleCell, 1.5, 5.5, 0.5, -70.42, 'midd', [5,10])
    pp_to_gc_x0 = [1000, 1000, 0.1]
    #loss(gc_to_in_x0, gc_to_in_args)
    res_pp_gc = minimize(loss, pp_to_gc_x0, pp_to_gc_args, method='Nelder-Mead', options={'maxfev': 3000, 'maxiter': 2000})
