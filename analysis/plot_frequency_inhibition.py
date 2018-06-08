# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import matplotlib.pyplot as plt

# Setup some parameters given by paradigm_frequency_inhibition.py
stim_delay = 100  # ms
dt = 0.01  # ms
stim_dtp = stim_delay / dt

data_path = "C:\\Users\\Daniel\\pyDentateData\\frequency_inhibition_data\\50Hz\\"
save_path = "C:\\Users\\DanielM\\Repos\\pyDentateData\\frequency_inhibition_data\\"
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f and '40' in f and not '50' in f]
hertz = 50
data = np.load(data_path + data_files[0])['arr_0']
for x in data_files[0:len(data_files)-1]:
    data = np.concatenate((data,np.load(data_path + x)['arr_0']))
    #data = np.load(data_path + x)['arr_0']
    #plt.plot(np.arange(0,500.02, 0.01), data.mean(axis=0))
    #for y in data:
    #    plt.plot(np.arange(0,500.02, 0.01), y)
