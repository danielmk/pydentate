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
data_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data_local_input_revised\\seed10000\\input_patterns_seed_10000\\"
save_path = data_path
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f and not 'norm' in f and not 'trifilt' in f]
data_files.sort()

data_files = data_files[0:25]

data_list = []
for x in data_files:
    data_list.append(np.load(data_path + x)['arr_0'])

row_idx_start = 0
row_idx_stop = 25
# 376
len_bins = 1000

corrs = []

for row_idx, x in enumerate(data_list):
    corrs.append(analysis_main.inner_pearsonr(x,len_bins))
        
np.savetxt(save_path + "1_inner-correlation_matrix_len-bin_" + str(len_bins) +  ".txt", corrs, delimiter="\t")