# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import analysis_main
from sklearn.preprocessing import normalize

# Setup some parameters given by paradigm_frequency_inhibition.py
stim_delay = 100  # ms
dt = 0.01  # ms
stim_dtp = stim_delay / dt

data_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data_local_input_revised\\scale1000\\net_tunedrev\\"
save_path = data_path
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.pydd' in f and not '.npz' in f]

for x in data_files:
    print("WHAT")
    data = shelve.open(data_path + x)
    curr_spike_data = analysis_main.time_stamps_to_signal(data['net_tunedrev.TunedNetwork']['populations'][0]['ap_time_stamps'], dt_signal=0.1, t_start=0, t_stop=600)
    np.savez(save_path + x + '_spike_data', curr_spike_data)
    curr_spike_data_conv = analysis_main.tri_filter(curr_spike_data, 200)
    np.savez(save_path + x + '_spike_data_convolved', curr_spike_data_conv)