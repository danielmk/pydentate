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

data_path = "C:\\Users\\Daniel\\pyDentateData\\tuning\\revised\\spatial_inhibition_data\\globalrev\\"
save_path = "C:\\Users\\Daniel\\pyDentateData\\tuning\\revised\\spatial_inhibition_data\\globalrev\\"
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.npz' in f and '2' in f[61]]
data_files = data_files[1::2]
#data_files=list(np.array(data_files)[[1,3,9,14,16,18,20,22,24,26]])
#data_files = data_files[0::2]
data = np.load(data_path + data_files[0])['arr_0']
stim_delay = 50
#intervals = np.arange(stim_delay/dt, (stim_delay/dt) + (interval/dt)*10, interval/dt)

for x in data_files[0:len(data_files)-1]:
    data = np.concatenate((data,np.load(data_path + x)['arr_0']))

peaks=[]
for x in data:
    x = x - x[int(30/dt):int(50/dt)].mean()
    data_peaks = []
    data_peaks.append(x[int(50/dt):int(100/dt)].max())
    #plt.plot(data_peaks, marker ='o')
    peaks.append(data_peaks)
    
peaks = np.array(peaks)
peaks[np.argwhere(peaks<0)]=0
peaks = peaks.flatten()
peaks = peaks.reshape((peaks.shape[0]/20,20))

np.savetxt(save_path + "paradigm-spatial-inhibition_peaks_n-cells_.txt", peaks, delimiter="\t")
