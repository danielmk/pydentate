# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 10:25:41 2018

@author: daniel
"""
import shelve
import numpy as np
import matplotlib.pyplot as plt

#Home PC
#directory = "C:\\Users\\daniel\\repos\\pyDentate\paradigm_pattern-separation_saves_2018-03-11\\"
#Office PC
directory = "C:\\Users\\DanielM\\Repos\\pyDentate\\"
file_name = "net_tuned.TunedNetwork_run_1"

shelve_files = []
active_cells = []
for run in range(1):
    shelve_files.append(shelve.open(directory + file_name + str(run)))
    ap_n_array = np.array(shelve_files[run]['net_tuned.TunedNetwork']['populations'][0]['ap_number'])
    active_cells.append(len(np.argwhere(ap_n_array)))

# Percent Active Cells
for x in shelve_files:
    n_active_array = np.array(x['net_tuned.TunedNetwork']['populations'][0]['ap_number'])
    n_active = len(np.argwhere(n_active_array))
    perc_active = (n_active / 2000.0) * 100
    print(perc_active)