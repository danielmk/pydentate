# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import matplotlib.pyplot as plt

sim_matrix_path = 'C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data\\net_global\\'
sim_matrix_fname = 'dot-product-norm_matrix_net_global_complete.txt'
input_matrix = np.loadtxt(sim_matrix_path + sim_matrix_fname, delimiter = '\t')

linear_similarity = np.array([])
for x in range(376):
    linear_similarity = np.concatenate((linear_similarity, input_matrix[x,x:input_matrix.shape[0]]))

np.savetxt(sim_matrix_path + 'linearized_similarity_matrix.txt', linear_similarity, delimiter = '\t')