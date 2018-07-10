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
from sklearn.preprocessing import normalize
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt = mpl.pyplot

save_dir = "C:\\Users\\Daniel\\pyDentateData\\input_pattern_plots\\"

seed = 10000

# Generate the temporal patternss
np.random.seed(seed)
temporal_patterns = inhom_poiss()

runs = range(25)

patterns = []
for x in runs:
    #curr_pattern = np.zeros(temporal_patterns.shape, dtype = np.ndarray)
    #curr_pattern = np.array([] * temporal_patterns.shape[0], dtype = np.ndarray)
    curr_pattern = np.empty((temporal_patterns.shape[0],), dtype=object)
    curr_pattern[:] = [[] * len(curr_pattern)] 
    curr_pattern[x:x+24] = temporal_patterns[x:x+24]
    if not bool(np.array(curr_pattern[0]).any()):
        curr_pattern[0] = [0]
    patterns.append(curr_pattern)

for idx, x in enumerate(patterns):
    plt.figure()
    plt.eventplot(x[0:50])
    plt.xlim((0,600))
    plt.savefig(save_dir + "input_pattern_eventplot_run" + str(idx).zfill(2) + ".pdf")
    plt.savefig(save_dir + "input_pattern_eventplot_run" + str(idx).zfill(2) + ".eps")
    plt.savefig(save_dir + "input_pattern_eventplot_run" + str(idx).zfill(2) + ".png")
    plt.close()