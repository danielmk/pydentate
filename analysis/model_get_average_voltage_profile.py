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
data_file = "net_tunedrev.TunedNetwork_run_scale_000_1000.pydd"
data = shelve.open(data_path + data_file)

avg_gcs = np.mean(data['net_tunedrev.TunedNetwork']['populations'][0]['v_records'], axis = 0)
avg_mcs = np.mean(data['net_tunedrev.TunedNetwork']['populations'][1]['v_records'], axis = 0)
avg_bcs = np.mean(data['net_tunedrev.TunedNetwork']['populations'][2]['v_records'], axis = 0)
avg_hcs = np.mean(data['net_tunedrev.TunedNetwork']['populations'][3]['v_records'], axis = 0)

np.savetxt(save_path + data_file + '_avg_volt_trace_GC.txt', avg_gcs, delimiter='\t', newline = '\n')
np.savetxt(save_path + data_file + '_avg_volt_trace_MC.txt', avg_mcs, delimiter='\t', newline = '\n')
np.savetxt(save_path + data_file + '_avg_volt_trace_BC.txt', avg_bcs, delimiter='\t', newline = '\n')
np.savetxt(save_path + data_file + '_avg_volt_trace_HC.txt', avg_hcs, delimiter='\t', newline = '\n')
