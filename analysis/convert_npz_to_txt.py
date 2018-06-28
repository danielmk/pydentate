# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve

# Setup some parameters given by paradigm_frequency_inhibition.py
stim_delay = 100  # ms
dt = 0.01  # ms
stim_dtp = stim_delay / dt

data_path = "C:\\Users\\Daniel\\pyDentateData\\tuning\\revised\\frequency_inhibition_data\\"
save_path = data_path
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f and not '.txt' in f]

for x in data_files:
    curr_arr = np.load(data_path + x)['arr_0']
    np.savetxt(save_path + '.'.join(x.split('.')[:-1]) + '.txt', curr_arr, delimiter='\t', newline= '\n')
