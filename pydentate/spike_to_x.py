# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 09:46:39 2023

@author: Daniel
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


def units_times_to_binary(units, times, n_units=784, dt=0.0001, total_time=0.1):
    """Convert spike times and their presynaptic unit to binary matrix"""
    binary = np.zeros((n_units, int(total_time / dt)), dtype=bool)

    times_to_idc = np.array(times / dt, dtype = int)

    binary[units, times_to_idc] = 1

    return binary


def units_times_to_jagged(units, times, n_units=784):
    """Convert spike times and their presynaptic unit to binary matrix"""
    output = [[] for x in range(n_units)]

    for i, u in enumerate(units):
        # pdb.set_trace()
        if u < n_units:
            output[u].append(times[i])

    return np.array([np.array(x) for x in output], dtype=object)

def units_times_to_rate(units, times, n_units=784, dt=0.0001, total_time=0.1):
    """Convert spike times and their presynaptic unit to binary matrix"""
    binary = units_times_to_binary(units, times, n_units, dt, total_time)
    
    rate = binary.sum(axis=1) / total_time

    return rate


def units_times_to_first_spike(units, times, n_units=784, dt=0.0001, total_time=0.1):
    binary = units_times_to_binary(units, times, n_units, dt, total_time)
    
    time_to_first_spike = np.array([np.argwhere(x)[0, 0] * dt if x.any() else 0.1 for x in binary])
    
    return time_to_first_spike


def units_times_to_average_time(units, times, n_units=784, dt=0.0001, total_time=0.1):
    binary = units_times_to_binary(units, times, n_units, dt, total_time)
    
    avg_time = np.array([(np.argwhere(x)[:, 0] * dt).mean() if x.any() else 0.1 for x in binary])
    
    return avg_time

def spikes_to_synaptic_current(units, times, n_units=784, dt=0.0001, total_time=0.1, tr=0.002, td=0.015):
    
    binary = units_times_to_binary(units, times, n_units=n_units, dt=dt, total_time=total_time)
    psc = np.zeros(n_units)
    hm = np.zeros(n_units)
    steps = range(int(total_time / dt))
    
    current_output = []

    # pdb.set_trace()

    for i in steps:
        psc = psc * np.exp(-dt / tr) + hm * dt
        
        # Integrate the current
        # hm = hm * np.exp(-dt / td) + JD * (len(index) > 0) / (tr * td)
        hm = hm * np.exp(-dt / td) + binary[:, i] / (tr * td)

        current_output.append(psc)
        
    return current_output

def spike_trains_to_units_times(spike_trains):
    units = []

    for idx, st in enumerate(spike_trains):
        units.append([idx] * len(st))
    units = np.concatenate(units)
    times = np.concatenate(spike_trains)
    # pdb.set_trace()
    idc_sorted = np.argsort(times)
    units = units[idc_sorted]
    times = times[idc_sorted]
    return units, times
        
        
        
        
        