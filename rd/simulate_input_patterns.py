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
from pyDentate.burst_generator_inhomogeneous_poisson import inhom_poiss
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# Generate the temporal patterns
seed = 10001
np.random.seed(seed)
temporal_patterns = inhom_poiss(rate=10)

runs = range(0,23,3)

temporal_patterns = analysis_main.time_stamps_to_signal(temporal_patterns, 0.1, 0, 600)
pattern_shape = temporal_patterns.shape

patterns = []
for x in runs:
    curr_pattern = np.zeros(temporal_patterns.shape)
    curr_pattern[x:x+24,:] = temporal_patterns[x:x+24,:]
    patterns.append(curr_pattern)

patterns = np.array(patterns)
patterns = patterns.sum(axis=2)
corrs = []
for idx_r, x in enumerate(patterns):
    for idx_c in range(1+idx_r, len(patterns)):
        curr_r = pearsonr(patterns[idx_r],patterns[idx_c])
        corrs.append(curr_r[0])
        
plt.figure()
plt.hist(corrs, bins=np.arange(0,1,0.1), range = (0, 1))
plt.xlabel("input corr")
plt.ylabel("# input patterns")
plt.title("Seed_ " + str(seed) + " Runs_ " + str(list(runs)))
plt.savefig("C:\\Users\\Daniel\\repos\\pyDentate\\rd\\" + "Seed_ " + str(seed) + " Runs_ " + str(list(runs)))

"""
np.savez(save_path + file_prefix + str(x).zfill(3), curr_pattern)
curr_pattern_trifilt[x:x+24,:] = temporal_patterns_trifilt[x:x+24,:]
np.savez(save_path + file_prefix + str(x).zfill(3) + '_trifilt', curr_pattern_trifilt)
curr_pattern_norm[x:x+24,:] = temporal_patterns_norm[x:x+24,:]
np.savez(save_path + file_prefix + str(x).zfill(3) + '_norm', curr_pattern_norm)
"""