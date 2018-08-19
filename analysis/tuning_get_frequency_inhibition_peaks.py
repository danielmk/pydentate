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

data_path = "C:\\Users\\Daniel\\pyDentateData\\tuning\\revised\\frequency_inhibition_data\\"
save_path = data_path
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f and '100' in f and 'HC' in f]
data_files = "paradigm-frequency-inhibition_interval_n-cells_1000_40_GCs_net_nonfacilitatingrev.txt"
interval = 1000 # in kHz because ms
data = np.loadtxt(data_path + data_files)
stim_delay = 100
intervals = np.arange(stim_delay/dt, (stim_delay/dt) + (interval/dt)*10, interval/dt)

for x in data_files[0:len(data_files)-1]:
    data = np.concatenate((data,np.load(data_path + x)['arr_0']))

peaks=[]
for x in data:
    x = x - x[int(80/dt):int(100/dt)].mean()
    data_peaks = []
    for start in intervals:
        print(start)
        data_peaks.append(x[int(start):int(start+(20/dt))].min())
    #plt.plot(data_peaks, marker ='o')
    peaks.append(data_peaks)
    
peaks = np.array(peaks)

norm_peaks = []
for idx, x in enumerate(peaks):
    x = (x/x[0])*100
    norm_peaks.append(x)
        #plt.plot(x, marker = 'o')
    
norm_peaks = np.array(norm_peaks)

np.savetxt(save_path + "paradigm-frequency-inhibition_norm-peaks_interval_n-cells_" + str(interval) + ".txt", norm_peaks, delimiter="\t")
np.savetxt(save_path + "paradigm-frequency-inhibition_peaks_interval_n-cells_" + str(interval) + ".txt", peaks, delimiter="\t")