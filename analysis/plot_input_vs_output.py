# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import matplotlib.pyplot as plt

# Setup some parameters given by paradigm_frequency_inhibition.py
stim_delay = 100  # ms
dt = 0.01  # ms
stim_dtp = stim_delay / dt

input_path = 'C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data\\input_patterns\\dot-product-norm_matrix_idx_row-start_0-376.txt'
output_path = 'C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data\\net_tuned\\dot-product-norm_matrix_complete_net_tuned_output.txt'

input_matrix = np.loadtxt(input_path, delimiter = '\t')
output_matrix = np.loadtxt(output_path, delimiter = '\t')

input_array = input_matrix.flatten()
output_array = output_matrix.flatten()