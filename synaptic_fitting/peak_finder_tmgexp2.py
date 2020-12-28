#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# PEP8
"""
Created on Thu Nov  7 14:48:40 2019

@author: barisckuru
"""

import numpy as np
import os
import time
start_time = time.time()

""" peakfinder() function takes the averaged data of multiple cells stimulated
w diff freqs as well as std data for each. It normalizes the data according
to first signal peak.Then it splits the data into ind pieces for each stimuli
and finds the peaks of each stimuli. It also finds the idcs for signal peaks
and uses it to extract stds in those points, also converts stds into sems.
Finally it returns normalized data, peaks in the data and sems of peaks. """


def peakfinder(data_path):
    load_1 = []
    peaks = []
    norm = []
    std = []
    peaks_std = []
    'LOAD THE EXPERIMENTAL DATA'
    data_path = data_path
    # Open and load the data in the same directory for diff freqs
    for files in os.walk(data_path):
        for file in files[2]:
            if '1hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                freq_1 = load_1['mean_arr']
            if '10hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                freq_10 = load_1['mean_arr']
            if '30hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                freq_30 = load_1['mean_arr']
            if '50hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                freq_50 = load_1['mean_arr']
    loads = [freq_1, freq_10, freq_30, freq_50]
    # Open and load the STD data in the same directory for diff freqs
    for files in os.walk(data_path):
        for file in files[2]:
            if '1hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                err_1 = load_1['std_arr']
            if '10hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                err_10 = load_1['std_arr']
            if '30hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                err_30 = load_1['std_arr']
            if '50hz' in file:
                curr_path = data_path + file
                load_1 = np.load(curr_path)
                err_50 = load_1['std_arr']
    errs = [err_1, err_10, err_30, err_50]

    'PEAK FINDER'
    # stimuli length is variable for diff freqs
    for i in range(len(loads)):
        # G data and std data for each freq were taken as "data" and "stds"
        data = loads[i]
        stds = errs[i]
        # stim_len is the range between stimuli
        # sig_len is only the length of the signal
        if loads[i] is freq_1:
            sig_len = 900
            stim_len = 20100
        elif loads[i] is freq_10:
            sig_len = 600
            stim_len = 2000
        elif loads[i] is freq_30:
            sig_len = 500
            stim_len = 666
        elif loads[i] is freq_50:
            sig_len = 300
            stim_len = 400
        stim_num = 10
        # Data is conductance, negative values converted to positive
        data = -data
        # first increase was found, can be caused either artifact or signal
        # artifact in these dataset and used to define staring point
        first_idc = np.argwhere(data > 50)
        # to eliminate artifact, first idc shifted 40 dps
        stim_start = first_idc[0] + 40
        stim_end = stim_start+stim_num*stim_len
        start_idcs = np.arange(stim_start, stim_end, stim_len)
        end_idcs = start_idcs + sig_len
        idcs = np.concatenate((start_idcs, end_idcs))
        idcs = np.sort(idcs)
        idc_len = len(idcs)
        # first sig isolated and signal peak was used for nowmalization
        first_sig = data[idcs[0]:idcs[1]]
        first_peak = max(first_sig)
        data = data/first_peak
        # intial time was cut
        data_cut = data[(idcs[0]-10000):(idcs[idc_len-1]+500)]
        # normalized & reshaped data appended
        norm.append(data_cut)
        # standart deviation was normalized
        stds = stds/first_peak
        stds_cut = stds[(idcs[0]-10000):(idcs[idc_len-1]+500)]
        std.append(stds_cut)
        # data and stds were splitted to take signals individually
        split_data = np.split(data, idcs)
        split_std = np.split(stds, idcs)
        split_data = split_data[1:len(split_data):2]
        split_std = split_std[1:len(split_std):2]
        split_data = np.array(split_data)
        split_std = np.array(split_std)
        # peaks and peak indices were found
        peaks_data = np.amax(split_data, axis=1)
        peaks.append(peaks_data)
        peaks_idc_split = np.argmax(split_data, axis=1)
        peaks_stds = []
        for i in range(len(peaks_idc_split)):
            peaks_stds.append((split_std[i][peaks_idc_split[i]]))
        peaks_std.append(peaks_stds)
    return peaks, norm, peaks_std
