# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve

# Setup some parameters given by paradigm_frequency_inhibition.py
stim_delay = 100  # ms
dt = 0.01  # ms
stim_dtp = stim_delay / dt

data_path = "C:\\Users\\DanielM\\Repos\\pyDentateData\\frequency_inhibition_data\\"
save_path = "C:\\Users\\DanielM\\Repos\\pyDentateData\\frequency_inhibition_data\\"
data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.pydd' in f]

for x in data_files:
    interval = int(x.split('_')[8].split('.')[0][1:3])
    data = shelve.open(data_path + x)
    split_name_current = x.split('.')
    split_name_peaks = list(split_name_current)
    split_name_current[1] = split_name_current[1] + '_current'
    split_name_peaks[1] = split_name_peaks[1] + '_peaks'
    name_current = '.'.join(split_name_current)
    name_peaks = '.'.join(split_name_peaks)
    np.savez(save_path + name_current, np.array(data[data.keys()[0]]['populations'][0]['VClamps_i']))

# Calculate spatial IPSC plot
bl_times = np.array([40, 50])  # in ms
IPSC_times = np.array([50, 90])  # in ms
bl_dtps = bl_times / sampling_period
IPSC_dtps = IPSC_times / sampling_period

IPSCs = []
for cell_i in nw.populations[0].VClamps_i:
    trace = cell_i.as_numpy()
    bl = trace[int(bl_dtps[0]):int(bl_dtps[1])].mean()
    peak_IPSC = trace[int(IPSC_dtps[0]):int(IPSC_dtps[1])].max()
    IPSCs.append(peak_IPSC - bl)
spatial_plot = plt.figure()
plt.plot(cells_to_measure, IPSCs)
plt.xlabel("Granule Cell #")
plt.ylabel("Peak IPSC (nA)")
full_file_path = save_dir + '\\' + 'run_' + str(run) + '_spatial_IPSC_plot_nw'
spatial_plot.savefig(full_file_path + ".pdf", dpi = 300, format ='pdf')
spatial_plot.savefig(full_file_path + ".eps", dpi = 300, format ='eps')