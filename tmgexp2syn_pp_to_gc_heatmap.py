# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 20:25:56 2019

@author: Daniel
"""

import numpy as np
from mossycell_cat import MossyCell
from basketcell import BasketCell
from granulecell import GranuleCell
import os
from neuron import h, gui
from scipy.optimize import minimize
import time
import pdb
from pyDentate.tmgexp2syn_fitter import loss, simulate

data_path = (
    "C:\\Users\\Daniel\\Dropbox\\02_MEC Project\\003_Antidromic "
    "electrically evoked EPSCs in hilar cells, voltage clamp\\"
)

pp_to_gc = np.load(data_path + "pp_to_gc_data_full.npz")
pp_to_gc_vec = np.concatenate(
    (pp_to_gc['peaks_norm_05'], pp_to_gc['peaks_norm_10'])
)
pp_to_gc_args = (pp_to_gc_vec, GranuleCell, 1.5, 5.5, 0.5, -70.42, 'midd', [5,10])
pp_to_gc_x0 = [10, 1000, 0.1]

u0 = 0.5
tau_facil_space = np.arange(0,1000,10)
tau_rec_space = np.arange(0,1000,10)
output_matrix = np.zeros((tau_facil_space.shape[0], tau_rec_space.shape[0]))
"""
for idx1, x in enumerate(tau_facil_space):
    for idx2, y in enumerate(tau_rec_space):
        sim = loss([x, y, u0], *pp_to_gc_args)
        output_matrix[idx1, idx2] = sim
    print(idx1)
"""

result1 = np.array([2.22590184e+13, 3.12395254e+02, 8.15226996e-01])
result1_sim = simulate(result1[0],
                       result1[1],
                       5,
                       pp_to_gc_args[1],
                       pp_to_gc_args[2],
                       pp_to_gc_args[3],
                       result1[2],
                       0.5,
                       pp_to_gc_args[5],
                       pp_to_gc_args[6])
