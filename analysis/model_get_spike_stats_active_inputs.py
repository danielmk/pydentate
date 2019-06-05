# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 10:25:41 2018

@author: daniel
"""
import shelve
import numpy as np
import matplotlib.pyplot as plt
import os
import pdb

#Home PC
#directory = "C:\\Users\\daniel\\repos\\pyDentate\paradigm_pattern-separation_saves_2018-03-11\\"
#Office PC
#directory = "Y:\\DanielM\\023_Dentate Gyrus Model\\paradigm_spatial-inhibition\\"
#Dropbox
def get_spike_stats(data_path):
    data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f and not 'norm' in f and not 'trifilt' in f]
    data_files.sort()

    # Get to BasketCell Connection
    n_spike_list = []
    for x in data_files:
        print(x)
        #curr_data = shelve.open(data_path + x)
        curr_data = np.load(data_path+x, mmap_mode='r')['arr_0']
        active_patterns = np.argwhere(np.any(curr_data, axis=1))[:,0]
        n_spikes = curr_data[active_patterns,:].sum(axis=1)
        n_spikes = np.append(n_spikes, n_spikes.mean())
        n_spike_list.append(n_spikes)

    np.savetxt(data_path + "1_n_spikes.txt", np.array(n_spike_list), delimiter = '\t')

if __name__ == '__main__':
    data_path = "Z:\\pyDentate\\pyDentateData\\pattern_separation_data_local_30Hz_input\\seed10006\\input_patterns\\"
    get_spike_stats(data_path)

