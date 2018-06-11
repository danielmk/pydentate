# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import os

# Setup some parameters given by paradigm_frequency_inhibition.py
stim_delay = 100  # ms
dt = 0.01  # ms
stim_dtp = stim_delay / dt

data_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data\\"
save_path = data_path
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.pydd' in f]

for x in data_files:
    split_name = x.split('_')
    num_split = x.split('_')[3].split('.')
    num_split[0] = num_split[0].zfill(3)
    split_name[3] = '.'.join(num_split)
    new_name = '_'.join(split_name)
    os.rename(data_path + x, data_path + new_name)