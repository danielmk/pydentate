# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import os
import analysis_main

# Setup some parameters given by paradigm_frequency_inhibition.py
stim_delay = 100  # ms
dt = 0.01  # ms
stim_dtp = stim_delay / dt

data_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data\\input_patterns\\"
save_path = data_path
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f]
data_files.sort()

corr_matrix = np.zeros((len(data_files), len(data_files)))

data_list = []
for x in data_files:
    data_list.append(np.load(data_path + x)['arr_0'])
    
for row_idx, x in enumerate(data_list):
    for col_idx, y in enumerate(data_list):
        corr_matrix[row_idx,col_idx]=analysis_main.correlate_signals(x,y)