# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 10:25:41 2018

@author: daniel
"""
import shelve
import numpy as np

directory = "C:\\Users\\daniel\\repos\\pyDentate\paradigm_pattern-separation_saves_2018-03-11\\"
file_name = "net_tuned.TunedNetwork_run_"

shelve_files = []
active_cells = []
for run in range(30):
    shelve_files.append(shelve.open(directory + file_name + str(run)))
    ap_n_array = np.array(shelve_files[run]['net_tuned.TunedNetwork']['populations'][0]['ap_number'])
    active_cells.append(len(np.argwhere(ap_n_array)))