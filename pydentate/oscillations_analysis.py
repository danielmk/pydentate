# -*- coding: utf-8 -*-
"""
This module contains the implementation of coherence and power measures used
throughout the study.

@author: danielmk
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import os
import sys
import uuid
import pdb
import shelve
from scipy.signal import windows, convolve
from scipy.ndimage import convolve
from scipy.spatial import distance_matrix
from numpy.fft import fft, fftshift, fftfreq

def pairwise_coherence(binary_spiketrain, bin_size=1):
    """This is the measure used in Wang & Buzsaki 1996.
    """
    
    binary_spiketrain = binary_spiketrain[:, :(binary_spiketrain.shape[1] // bin_size) * bin_size]
    
    binned_spiketrain = binary_spiketrain.reshape(binary_spiketrain.shape[0], binary_spiketrain.shape[1] // bin_size, bin_size).sum(axis=-1)
    
    cross_sum = np.matmul(binned_spiketrain, binned_spiketrain.T)
    
    summed_output = (binned_spiketrain).sum(axis=1) 
    
    cross_output = np.sqrt(np.matmul(summed_output[:, None], summed_output[None, :]))
    
    coherence_matrix = cross_sum / cross_output

    ut_idc = np.triu_indices(binary_spiketrain.shape[0], 1)

    return coherence_matrix[ut_idc]

def linearity_analysis(binary_spiketrain, dt, duration=2, n_points=30):
    """Estimate the area between a line that goes through the points
    (0, 0) and (tau_(1/f_u), k(tau_(1/f_u))) and the function k(tau), where
    k is the pairwise_coherence.
    """
    lp0 = (0, 0)

    n_spikes = binary_spiketrain.sum(axis=1)
    n_spikes_firing = n_spikes[n_spikes > 0]
    avg_freq = n_spikes_firing.mean() / duration  # Convert duration to ms
    peak_tau_ms = 1 / (avg_freq / 1000)
    
    peak_bin_size = int(peak_tau_ms / dt)
    
    peak_tau_ms_coherence = np.nanmean(pairwise_coherence(binary_spiketrain, peak_bin_size))

    line_slope = peak_tau_ms_coherence / peak_tau_ms
    
    x_values_ms = np.linspace(0, peak_tau_ms, n_points)
    
    y_values_line = x_values_ms * line_slope

    y_values_k = np.array([np.nanmean(pairwise_coherence(binary_spiketrain, int(x / dt))) if x != 0 else 0 for x in x_values_ms])
    
    area_estimate = np.array([(x_values_ms[i] - x_values_ms[i - 1]) * (y_values_k[i - 1] - y_values_line[i]) +
                              (1/2) * ((x_values_ms[i] - x_values_ms[i - 1]) * (y_values_line[i] - y_values_line[i - 1])) +
                              (1/2) * ((x_values_ms[i] - x_values_ms[i - 1]) * np.abs(y_values_k[i] - y_values_k[i - 1]))
                              for i in range(1, n_points)])
                              
    return area_estimate.sum(), area_estimate.sum() / peak_tau_ms

def shuffling_coherence(times, dt, duration, n_shuffles=100):
    """Coherence through 

    """
    n_timepoints = int((duration * 1000) / dt)
    
    spike_counts = np.unique(times, return_counts=True)
    
    spike_occurences_indices = np.array(spike_counts[0] / dt, dtype=int)
    
    spike_counts_full = np.zeros(n_timepoints+1)
    
    spike_counts_full[spike_occurences_indices] = spike_counts[1]
    
    gaussian_sigma = 1 / dt
    gaussian_width = gaussian_sigma * 10
    gaussian_kernel = windows.gaussian(gaussian_width, gaussian_sigma)
    gaussian_kernel = gaussian_kernel / gaussian_kernel.sum()
    
    spike_counts_full_convolved = convolve(spike_counts_full, gaussian_kernel)
    
    # Shuffle spikes
    random_times = np.array([np.random.uniform(0, duration * 1000, times.shape[0]) for x in range(n_shuffles)])
    
    bin_edges = np.arange(0, duration * 1000 + 0.2, dt)
    
    random_spike_counts_full = np.array([((random_times > bin_edges[i - 1]) & (random_times < bin_edges[i])).sum() for i in range(1, bin_edges.shape[0])], dtype=float) / n_shuffles
    # pdb.set_trace()
    random_spike_counts_full_convolved = convolve(random_spike_counts_full, gaussian_kernel)

    coherence = np.abs((spike_counts_full_convolved - random_spike_counts_full_convolved)).sum()

    return spike_counts_full_convolved, random_spike_counts_full_convolved

"""
def k_(binary_spiketrain, bin_size=1):
    
    binary_spiketrain = binary_spiketrain[:, :(binary_spiketrain.shape[1] // bin_size) * bin_size]
    
    binned_spiketrain = binary_spiketrain.reshape(binary_spiketrain.shape[0], binary_spiketrain.shape[1] // bin_size, bin_size).sum(axis=-1)
    
    cross_sum = np.matmul(binned_spiketrain, binned_spiketrain.T)
    
    summed_output = (binned_spiketrain ** 2).sum(axis=1) 
    
    cross_output = np.sqrt(np.matmul(summed_output[:, None], summed_output[None, :]))
    
    coherence_matrix = cross_sum / cross_output

    ut_idc = np.triu_indices(binary_spiketrain.shape[0], 1)
    # pdb.set_trace()
    return coherence_matrix[ut_idc]
"""
