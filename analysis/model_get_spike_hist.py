# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import matplotlib.pyplot as plt

data_path = "C:\\Users\\Daniel\\pyDentateData\\example\\"
save_path = data_path
data_file = "net_tunedrev.TunedNetwork_data_paradigm_local-pattern-separation-30Hz_run_scale_seed_000_1000_10000.pydd"
data = shelve.open(data_path + data_file)

spikes = data['net_tunedrev.TunedNetwork']['populations'][0]['ap_time_stamps']
spikes_flat = np.concatenate(spikes)
myhist = plt.hist(spikes_flat, bins = np.arange(0,600,5))