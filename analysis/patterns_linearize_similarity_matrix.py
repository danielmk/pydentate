# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import matplotlib.pyplot as plt

sim_matrix_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data_local_repeated_input_revised\\seed10000\\scale1000\\net_nonfacilitatingrev\\"
sim_matrix_fname = '1_leutgeb-measure_matrix_len-bin_6000.txt'
input_matrix = np.loadtxt(sim_matrix_path + sim_matrix_fname, delimiter = '\t')

linear_similarity = np.array([])
for x in range(25):
    linear_similarity = np.concatenate((linear_similarity, input_matrix[x,x:input_matrix.shape[0]]))

np.savetxt(sim_matrix_path + '1_linearized_leutgeb-measure_matrix_len-bin_6000.txt', linear_similarity, delimiter = '\t')