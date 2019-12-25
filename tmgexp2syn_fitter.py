# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 21:57:24 2019

@author: Daniel
"""

import numpy as np
from pyDentate.mossycell_cat import MossyCell
from basketcell import BasketCell
import os
from neuron import h, gui
from scipy.optimize import minimize
import time

# from peak_finder_tmgexp2 import peakfinder
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
    peaks_sem : ndarray
        Standard error of means of the peaks.
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
    peaks_sem = peaks_sem / peaks[0]

    return peaks_norm, peaks_sem


def loss(x, *args):
    """Loss function for use with scipy.optimize.minimize. x is a 1d-array
    that presents the variables that are to be optimized and *args contains all
    other parameters needed to specify the function.
    
    Parameters
    ----------
    x : 1d-array
        x[0], tau_facil. x[1], tau_rec.j x[2], u0.
    *args : tuple
        args[0], vectorized peaks from the data in order of frequencies
        args[1], celltype.
        args[2], tau_1.
        args[3], tau_2.
        args[4], sampling_interval.
        args[5], v_clamp.
        args[6], sec.
    """
    # Run the simuations
    sim_01Hz = simulate(
        x[0],
        x[1],
        1,
        args[1],
        tau_1=args[2],
        tau_2=args[3],
        u0=x[2],
        sampling_interval=args[4],
        v_clamp=args[5],
        sec=args[6]
    )

    sim_10Hz = simulate(
        x[0],
        x[1],
        10,
        args[1],
        tau_1=args[2],
        tau_2=args[3],
        u0=x[2],
        sampling_interval=args[4],
        v_clamp=args[5],
        sec=args[6]
    )

    sim_30Hz = simulate(
        x[0],
        x[1],
        30,
        args[1],
        tau_1=args[2],
        tau_2=args[3],
        u0=x[2],
        sampling_interval=args[4],
        v_clamp=args[5],
        sec=args[6]
    )

    sim_50Hz = simulate(
        x[0],
        x[1],
        50,
        args[1],
        tau_1=args[2],
        tau_2=args[3],
        u0=x[2],
        sampling_interval=args[4],
        v_clamp=args[5],
        sec=args[6]
    )

    sim_peaks = np.concatenate((sim_01Hz[0], sim_10Hz[0], sim_30Hz[0], sim_50Hz[0]))

    mean_squared_error = np.square(sim_peaks - args[0]).mean()

    return mean_squared_error


if __name__ == "__main__":
    data_path = (
        "C:\\Users\\Daniel\\Dropbox\\02_MEC Project\\003_Antidromic "
        "electrically evoked EPSCs in hilar cells, voltage clamp\\"
    )
    gc_to_in = np.load(data_path + "gc_to_in_avg_waveform_n-10.npz")
    gc_to_mc = np.load(data_path + "gc_to_mc_avg_waveform_n-10.npz")

    gc_to_in_01Hz = peakfinder(
        gc_to_in["mean_arr_01"],
        gc_to_in["sem_arr_01"],
        gc_to_in["stim_01"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    gc_to_in_10Hz = peakfinder(
        gc_to_in["mean_arr_10"],
        gc_to_in["sem_arr_10"],
        gc_to_in["stim_10"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    gc_to_in_30Hz = peakfinder(
        gc_to_in["mean_arr_30"],
        gc_to_in["sem_arr_30"],
        gc_to_in["stim_30"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    gc_to_in_50Hz = peakfinder(
        gc_to_in["mean_arr_50"],
        gc_to_in["sem_arr_50"],
        gc_to_in["stim_50"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    gc_to_in_vec = np.concatenate(
        (gc_to_in_01Hz[0], gc_to_in_10Hz[0], gc_to_in_30Hz[0], gc_to_in_50Hz[0])
    )
    gc_to_in_args = (gc_to_in_vec, BasketCell, 0.3, 0.6, 0.5, -70.42, 'proxd')
    gc_to_in_x0 = [3000, 10, 0.1]
    #loss(gc_to_in_x0, gc_to_in_args)
    res_gc_in = minimize(loss, gc_to_in_x0, gc_to_in_args, method='Nelder-Mead')

    
    gc_to_mc_01Hz = peakfinder(
        gc_to_mc["mean_arr_01"],
        gc_to_mc["sem_arr_01"],
        gc_to_mc["stim_01"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    gc_to_mc_10Hz = peakfinder(
        gc_to_mc["mean_arr_10"],
        gc_to_mc["sem_arr_10"],
        gc_to_mc["stim_10"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    gc_to_mc_30Hz = peakfinder(
        gc_to_mc["mean_arr_30"],
        gc_to_mc["sem_arr_30"],
        gc_to_mc["stim_30"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    gc_to_mc_50Hz = peakfinder(
        gc_to_mc["mean_arr_50"],
        gc_to_mc["sem_arr_50"],
        gc_to_mc["stim_50"],
        peak_win=(0.002, 0.016),
        peaks="negative",
        sampling_rate=20000,
    )
    
    gc_to_mc_vec = np.concatenate(
        (gc_to_mc_01Hz[0], gc_to_mc_10Hz[0], gc_to_mc_30Hz[0], gc_to_mc_50Hz[0])
    )
    
    # tau_1=0.5, tau_2=6.2 values gc to mc from the Santhakumar 2005 paper
    gc_to_mc_args = (gc_to_mc_vec, MossyCell, 0.5, 6.2, 0.5, -70.42, 'proxd')
    gc_to_mc_x0 = [3000, 25, 0.1]
    #loss(gc_to_in_x0, gc_to_in_args)
    res_gc_mc = minimize(loss, gc_to_mc_x0, gc_to_mc_args, method='Nelder-Mead')

    """
    # Data path msut be given here either for gcmc or gcin, or another dataset
    # unused one was commented out
    #data_path_gcmc = "/home/can/Downloads/gc_to_mc"
    data_path_gcmc = "C:\\Users\\Daniel\\Dropbox\\02_MEC Project\\BarisProject\\data"
    #data_path_gcin = "/home/can/Downloads/gc_to_in"
    peaks = peakfinder(data_path_gcmc)[0]
    #peaks = peakfinder(data_path_gcin)[0]
    taus = []
    'LOSS function'
    # with Mean Square Error for each freq and then avarage of all'
    
    
    def loss(x):
        tau_facil, tau_rec, u0 = x
        # taus and u0 are collected for each iteration of loss func
        taus.append(x)
        output1hz = simulate(x[0], x[1], 1, x[2])[0]
        Hz1 = peaks[0]
        output10hz = simulate(x[0], x[1], 10, x[2])[0]
        Hz10 = peaks[1]
        output30hz = simulate(x[0], x[1], 30, x[2])[0]
        Hz30 = peaks[2]
        output50hz = simulate(x[0], x[1], 50, x[2])[0]
        Hz50 = peaks[3]
    
        mse1 = (np.square(output1hz - Hz1)).mean(axis=None)
        mse10 = (np.square(output10hz - Hz10)).mean(axis=None)
        mse30 = (np.square(output30hz - Hz30)).mean(axis=None)
        mse50 = (np.square(output50hz - Hz50)).mean(axis=None)
        mse = (mse1 + mse10 + mse30 + mse50)/4
    
        return mse
    
    
    'OPTIMIZATION'
    
    # initial array to optimize is given here as x0
    x0 = np.array([500, 0, 0.1])
    # resulting MSE of loss function is minimized by Nelder-Mead Algorithm
    res = minimize(loss, x0, method='Nelder-Mead')
    # results are saved
    np.savez('gc_to_mc_opt', loss = res)
    # np.savez('gc_to_in_opt', loss = res)
    # time that optimizaation takes 
    end = time.time()
    opt_time = end-begin
    print (opt_time)
    """
