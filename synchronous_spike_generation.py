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

if __name__=='__main__':
    max_time = 2.4
    n_trial = 192
    rate = 10
    n_stim = 20
    offset = 0.1

    indices_unsorted = np.arange(0,n_trial)

    n_active = 150 + np.random.randint(0,42+1, n_stim)

    units_times = [[np.sort(np.random.permutation(indices_unsorted)[:n]), np.array([offset+idx*(1/rate)]*n)] for idx, n in enumerate(n_active)]

    units_arr = np.concatenate([x[0] for x in units_times], dtype=int)
    
    times_arr = np.concatenate([x[1] for x in units_times])

    shifts = np.arange(0.0, 0.04, 0.002)

    units_times_arr = [(units_arr, times_arr + np.random.normal(scale=x, size=times_arr.shape[0])) for x in shifts]

    pl_binary = [spike_to_x.units_times_to_binary(x[0], x[1]*1000, n_units=n_trial, dt=0.1, total_time=max_time*1000) for x in units_times_arr]

    linearity_measures = np.array([oscillations_analysis.linearity_analysis(x, dt=0.1, duration=max_time, n_points=30) for x in pl_binary])

    def exp_decay(t, tau, V):
        """Exponential decay helper function"""
        return V * np.e ** (-t / tau)

    # Create the synaptic kernel
    decay_tau = 0.0018
    t = np.arange(0, 0.020, 0.0001)
    kernel = exp_decay(t, decay_tau, 1)
    kernel = kernel / kernel.sum()

    binary_convolved = np.array([[np.convolve(x, kernel, 'full') for x in curr_binary ] for curr_binary in pl_binary])
    pr = np.array([np.corrcoef(x=curr_binary_convolved[:,:]) for curr_binary_convolved in binary_convolved])
    ut_idc = np.triu_indices(pr.shape[1], 1)
    mean_pearson_r = np.array([np.nanmean(curr_pr[ut_idc]) for curr_pr in pr])

    coherence = []
    for curr_binary in pl_binary:
        mean_frequency = (curr_binary.sum(axis=1) / max_time).mean()
        # print(mean_frequency)
        delta_t_1 = (0.1 / mean_frequency) * 1000
        curr_dt=0.1
        curr_coherence = oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(delta_t_1 / curr_dt)).mean()
        coherence.append(curr_coherence)

    coherence = np.array(coherence)

    plt.figure()
    plt.plot(coherence[::-1] / coherence.max(), marker='o')
    plt.plot(mean_pearson_r[::-1] / mean_pearson_r.max(), marker='o')
    plt.plot(linearity_measures[::-1,0] / linearity_measures[::-1,0].max(), marker='o')
    
    fig, ax = plt.subplots(3,1)
    ax[0].plot(coherence[::-1], marker='o')
    ax[1].plot(mean_pearson_r[::-1], marker='o')
    ax[2].plot(linearity_measures[::-1,0], marker='o')

    """
    shifts = np.arange(0.0,100,5)
    
    spike_trains = [cpp(10*pq.Hz,y,max_time * pq.s,shift=x) for x in shifts]
    
    ccp_sts = [[x.times for x in st] for st in spike_trains]

    ccp_units_times = [spike_to_x.spike_trains_to_units_times(x) for x in ccp_sts]
    
    pl_binary = [spike_to_x.units_times_to_binary(x[0], x[1]*1000, n_units=n_trial, dt=0.1, total_time=max_time*1000) for x in ccp_units_times]
    
    linearity_measures = [oscillations_analysis.linearity_analysis(x, dt=0.1, duration=max_time, n_points=30) for x in pl_binary]
    """

    
    
    
    # linearity_measure = [oscillations_analysis.linearity_analysis(binary_spiketrain, dt, duration=max_time, n_points=30)
    
    
    #oscillations_analysis.linearity_analysis(binary_spiketrain, dt, duration=2, n_points=30):
    #oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(bs/dt))