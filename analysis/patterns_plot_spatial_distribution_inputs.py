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
from scipy import stats
import matplotlib.pyplot as plt

seed = 10000
input_scale = 1000

# Generate the temporal patternss
np.random.seed(seed)
# Generate a gaussian probability density function
gauss_gc = stats.norm(loc=1000, scale=input_scale)
gauss_bc = stats.norm(loc=12, scale=(input_scale/2000.0)*24)
pdf_gc = gauss_gc.pdf(np.arange(2000))
pdf_gc = pdf_gc/pdf_gc.sum()
pdf_bc = gauss_bc.pdf(np.arange(24))
pdf_bc = pdf_bc/pdf_bc.sum()
# We hold the pdf constant. To randomize the centroid we reslice the GC indices
GC_indices = np.arange(2000)
start_idc_gc = np.random.randint(0, 1999, size=400)

PP_to_GCs = []
for x in start_idc_gc:
    curr_idc = np.concatenate((GC_indices[x:2000], GC_indices[0:x]))
    PP_to_GCs.append(np.random.choice(curr_idc, size=100, replace=False, p=pdf_gc))

PP_to_GCs = np.array(PP_to_GCs)
# Generate the PP -> BC mapping as above
BC_indices = np.arange(24)
start_idc_bc = np.array(((start_idc_gc/2000.0)*24), dtype=int)

PP_to_BCs = []
for x in start_idc_bc:
    curr_idc = np.concatenate((BC_indices[x:24], BC_indices[0:x]))
    PP_to_BCs.append(np.random.choice(curr_idc, size=1, replace=False, p=pdf_bc))

PP_to_BCs = np.array(PP_to_BCs)

aligned_targets = []
for x in PP_to_GCs:
    curr_PP = x.copy()
    curr_PP = (curr_PP - np.round(np.median(x))) + 1000
    curr_PP[curr_PP<0] = curr_PP[curr_PP<0] + 2000
    curr_PP[curr_PP>2000] = curr_PP[curr_PP>2000] - 2000
    aligned_targets.append(curr_PP)
    print(np.round(np.median(x)), curr_PP)
aligned_targets = np.array(aligned_targets)
plt.hist(aligned_targets.flatten(), bins = np.arange(-0.5,2000, 100))
plt.xlabel("# Granule Cell")
plt.ylabel("# Synapses received by GC")
plt.xlim((-0.5,2000))
plt.title("Seed: " + str(seed) + "; Scale: " + str(input_scale) + "; bin size: 100")


