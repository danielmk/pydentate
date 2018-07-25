# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 15:37:20 2018

@author: DanielM
"""

import numpy as np
import matplotlib.pyplot as plt
from elephant import spike_train_generation
from neo.core import AnalogSignal
import quantities as pq

def inhom_poiss_30Hz():
    sampling_interval = 0.0001 * pq.s
    max_rate = 100

    t = np.arange(0, 0.5, sampling_interval.magnitude)

    rate_profile = (np.sin(t*30*np.pi*2-np.pi/2) + 1) * max_rate / 2

    rate_profile_as_asig = AnalogSignal(rate_profile,
                                        units=1*pq.Hz,
                                        t_start=0*pq.s,
                                        t_stop=0.5*pq.s,
                                        sampling_period=sampling_interval)

    spike_trains = []
    for x in range(400):
        curr_train = spike_train_generation.inhomogeneous_poisson_process(rate_profile_as_asig)
        # We have to make sure that there is sufficient space between spikes.
        # If there is not, we move the next spike by 0.1ms
        spike_trains.append(curr_train)

    array_like = np.array([np.around(np.array(x.times)*1000, decimals=1) for x in spike_trains])
    for arr_idx in range(array_like.shape[0]):
        bad_idc = np.argwhere(np.diff(array_like[arr_idx]) == 0).flatten()
        bad_idc = bad_idc+1
        while bad_idc.any():
            for bad_idx in bad_idc:
                array_like[arr_idx][bad_idx] = array_like[arr_idx][bad_idx] + 0.1
            bad_idc = np.argwhere(np.diff(array_like[arr_idx]) == 0).flatten()
            bad_idc = bad_idc + 1

    return array_like

def inhom_poiss(rate=10, max_rate=100):
    sampling_interval = 0.0001 * pq.s

    t = np.arange(0, 0.5, sampling_interval.magnitude)

    rate_profile = (np.sin(t*rate*np.pi*2-np.pi/2) + 1) * max_rate / 2

    rate_profile_as_asig = AnalogSignal(rate_profile,
                                        units=1*pq.Hz,
                                        t_start=0*pq.s,
                                        t_stop=0.5*pq.s,
                                        sampling_period=sampling_interval)

    spike_trains = []
    for x in range(400):
        curr_train = spike_train_generation.inhomogeneous_poisson_process(rate_profile_as_asig)
        # We have to make sure that there is sufficient space between spikes.
        # If there is not, we move the next spike by 0.1ms
        spike_trains.append(curr_train)

    array_like = np.array([np.around(np.array(x.times)*1000, decimals=1) for x in spike_trains])
    for arr_idx in range(array_like.shape[0]):
        bad_idc = np.argwhere(np.diff(array_like[arr_idx]) == 0).flatten()
        bad_idc = bad_idc+1
        while bad_idc.any():
            for bad_idx in bad_idc:
                array_like[arr_idx][bad_idx] = array_like[arr_idx][bad_idx] + 0.1
            bad_idc = np.argwhere(np.diff(array_like[arr_idx]) == 0).flatten()
            bad_idc = bad_idc + 1

    return array_like

def hom_poiss(rate):
    rate = rate * pq.Hz

    spike_trains = []
    for x in range(400):
        curr_train = spike_train_generation.homogeneous_poisson_process(rate, 0.01*pq.s, 0.51 *pq.s)
        # We have to make sure that there is sufficient space between spikes.
        # If there is not, we move the next spike by 0.1ms
        spike_trains.append(curr_train)

    array_like = np.array([np.around(np.array(x.times)*1000, decimals=1) for x in spike_trains])
    for arr_idx in range(array_like.shape[0]):
        bad_idc = np.argwhere(np.diff(array_like[arr_idx]) < 0.1).flatten()
        bad_idc = bad_idc+1
        while bad_idc.any():
            for bad_idx in bad_idc:
                array_like[arr_idx][bad_idx] = array_like[arr_idx][bad_idx] + 0.1
            bad_idc = np.argwhere(np.diff(array_like[arr_idx]) < 0.1).flatten()
            bad_idc = bad_idc + 1

    return array_like