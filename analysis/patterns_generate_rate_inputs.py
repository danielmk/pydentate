# -*- coding: utf-8 -*-
"""
This script generates the input patterns as they are generated in paradigm_pattern_separation
and saves them to individual files for later processing.

@author: DanielM
"""

import os
import numpy as np
import shelve
import analysis_main
from pyDentate.burst_generator_inhomogeneous_poisson import inhom_poiss, hom_poiss
from sklearn.preprocessing import normalize

# Generate the temporal patterns
np.random.seed(0)
temporal_patterns = hom_poiss(100)

runs = range(25)
save_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data_rate\\input_patterns_100\\"
file_prefix = "input_patterns_rate_seed_100_0"

temporal_patterns = analysis_main.time_stamps_to_signal(temporal_patterns, 0.1, 0, 600)
temporal_patterns_trifilt = analysis_main.tri_filter(temporal_patterns, 200)
temporal_patterns_norm = normalize(temporal_patterns_trifilt, axis=1)
pattern_shape = temporal_patterns.shape

for x in runs:
    curr_pattern = np.zeros(temporal_patterns.shape)
    curr_pattern_norm = np.zeros(temporal_patterns.shape)
    curr_pattern_trifilt = np.zeros(temporal_patterns.shape)
    curr_pattern[x:x+24,:] = temporal_patterns[x:x+24,:]
    np.savez(save_path + file_prefix + str(x).zfill(3), curr_pattern)
    curr_pattern_trifilt[x:x+24,:] = temporal_patterns_trifilt[x:x+24,:]
    np.savez(save_path + file_prefix + str(x).zfill(3) + '_trifilt', curr_pattern_trifilt)
    curr_pattern_norm[x:x+24,:] = temporal_patterns_norm[x:x+24,:]
    np.savez(save_path + file_prefix + str(x).zfill(3) + '_norm', curr_pattern_norm)
