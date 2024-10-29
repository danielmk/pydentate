# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 09:13:44 2024

@author: Daniel
"""

import numpy as np
from pydentate import net_basket_cell_ring, neuron_tools, oscillations_analysis, spike_to_x
import matplotlib.pyplot as plt
from elephant.spike_train_generation import cpp
import scipy.stats
import quantities as pq

def temporally_correlated(max_time, n_trial, rate, regularity):
    """
    Generates a set of spike times with correlations within trains.

    Parameters:
    max_time (float): Maximum time for generating spike times.
    n_trial (int): Number of trials.
    rate (float): Average spike rate.
    regularity (float): Regularity parameter for gamma distribution.

    Returns:
    input_spiketimes (list of np.ndarray): List containing spike times for each trial.
    """
    
    input_spiketimes = []
    
    for trial_ind in range(n_trial):
        these_spiketimes = []
        running_time = np.random.exponential(1 / rate)
        n_spike = 0
        
        while running_time < max_time:
            n_spike += 1
            these_spiketimes.append(running_time)
            running_time += np.random.gamma(regularity, 1 / (rate * regularity))
        
        input_spiketimes.append(np.array(these_spiketimes))
    
    return input_spiketimes

def phase_locked(max_time, n_trial, rate, amplitude, num_phase):
    """
    Generates a set of spike times with sinusoidally varying rates.
    
    Parameters:
        max_time (float): Maximum time for simulation.
        n_trial (int): Number of trials.
        rate (float): Base firing rate.
        amplitude (float): Amplitude of the sinusoidal modulation.
        num_phase (int): Number of phases in the sinusoidal modulation.
    
    Returns:
        input_spiketimes (list of lists): List of spike times for each trial.
    """
    dt = 0.001
    input_spiketimes = []

    for trial_ind in range(n_trial):
        these_spiketimes = []
        n_spike = 0

        for time_ind in range(int(np.ceil(max_time / dt))):
            t_here = (time_ind + 0.5) * dt
            rate_here = rate * (1 + amplitude * np.sin(2 * np.pi * num_phase * t_here / max_time))
            p = np.random.rand()

            if p < (dt * rate_here):
                n_spike += 1
                these_spiketimes.append(t_here)

        input_spiketimes.append(np.array(these_spiketimes))

    return np.array(input_spiketimes, dtype=object)

if __name__=='__main__':
    max_time = 2
    n_trial = 192
    rate = 10
    regularity = 0.1
    amplitude = 1
    num_phase = 10
    # tc_st = temporally_correlated(max_time, n_trial, rate, regularity)
    # pl_st = phase_locked(max_time, n_trial, rate, amplitude, num_phase)
    """
    amplitudes = np.arange(0.1,50,1)

    pl_sts = [phase_locked(max_time, n_trial, rate, a, num_phase) for a in amplitudes]

    pl_units_times = [spike_to_x.spike_trains_to_units_times(x) for x in pl_sts]
    
    pl_binary = [spike_to_x.units_times_to_binary(x[0], x[1]*1000, n_units=n_trial, dt=0.1, total_time=max_time*1000) for x in pl_units_times]
    
    #linearity_measures = [oscillations_analysis.linearity_analysis(x, dt=0.1, duration=max_time, n_points=30) for x in pl_binary]
    
    units, times = spike_to_x.spike_trains_to_units_times(pl_sts[0])
    
    fig, ax = plt.subplots(2, 1)
    ax[0].eventplot(pl_sts[0])
    ax[1].eventplot(pl_sts[-1])
    """
    x=np.arange(1, n_trial+1)
    y = scipy.stats.norm.pdf(x, loc=n_trial/2, scale=1)
    y = y / y.sum()
    
    shifts = np.arange(0.0,100,5) * pq.ms
    
    spike_trains = [cpp(10*pq.Hz,y,max_time * pq.s,shift=x) for x in shifts]
    
    ccp_sts = [[x.times for x in st] for st in spike_trains]

    ccp_units_times = [spike_to_x.spike_trains_to_units_times(x) for x in ccp_sts]
    
    pl_binary = [spike_to_x.units_times_to_binary(x[0], x[1]*1000, n_units=n_trial, dt=0.1, total_time=max_time*1000) for x in ccp_units_times]
    
    linearity_measures = [oscillations_analysis.linearity_analysis(x, dt=0.1, duration=max_time, n_points=30) for x in pl_binary]
    
    plt.figure()
    plt.eventplot(ccp_sts[-1])
    
    
    
    # linearity_measure = [oscillations_analysis.linearity_analysis(binary_spiketrain, dt, duration=max_time, n_points=30)
    
    
    #oscillations_analysis.linearity_analysis(binary_spiketrain, dt, duration=2, n_points=30):
    #oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(bs/dt))