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
data_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data_local_input_revised\\seed10000\\scale1000\\net_globalrev\\"
save_path = data_path
data_files = [
    f
    for f in os.listdir(data_path)
    if os.path.isfile(os.path.join(data_path, f))
    and '.npz' in f
    and 'spike_data' in f
    and 'convolved' not in f
]

data_files.sort()

data_files = data_files[0:25]

corr_matrix = np.empty((len(data_files), len(data_files)))

data_list = []
for x in data_files:
    data_list.append(np.load(save_path + x)['arr_0'])
row_idx_start = 0
row_idx_stop = 25
# 376
len_bins = 6000

for row_idx, x in enumerate(data_list[row_idx_start:row_idx_stop]):
    for col_idx, y in enumerate(data_list[row_idx+row_idx_start:len(data_list)]):
        corr_matrix[row_idx+row_idx_start,col_idx+row_idx+row_idx_start]=analysis_main.sqrt_diff(x,y,len_bins)

np.savetxt(save_path + "1_sqrt-diff_matrix_len-bin_" + str(len_bins) +  ".txt", corr_matrix, delimiter="\t")