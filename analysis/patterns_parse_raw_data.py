# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 15:55:10 2019

@author: Daniel
"""

from scipy.stats import pearsonr
import numpy as np
import os 

directory = "C:\\Users\\Daniel\\ps_data\\"

all_files = os.listdir(directory)
all_files_split = np.array([np.array(x.split('_')) for x in all_files if os.path.isfile(directory + x)])
input_files = [x for x in all_files_split if x[0] == 'input']
output_files = [x for x in all_files_split if x[0] == 'output']
input_files_parameters = np.array([np.array(x[4:]) for x in all_files_split if x[0] == 'input'])
parameter_set = np.unique(input_files_parameters,axis=0)
files_parameter_one = np.array([x for x in all_files_split if np.array_equal(x[4:], parameter_set[0])])
for x in parameter_set:
    pass