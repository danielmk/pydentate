# -*- coding: utf-8 -*-
"""
Created on Sun May 20 11:35:03 2018

@author: DanielM
"""

import os
import numpy as np
import shelve
import analysis_main
from sklearn.preprocessing import normalize
from tempfile import TemporaryFile

def get_binary_traces(data_path):
    save_path = data_path
    data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.pydd' in f and not '.npz' in f]

    for x in data_files:
        if not os.path.isfile(data_path + '\\' + x +  '_spike_data.npz'):
            print(x)
            data = shelve.open(data_path + '\\' + x)
            curr_spike_data = analysis_main.time_stamps_to_signal(data[data.keys()[0]]['populations'][0]['ap_time_stamps'], dt_signal=0.1, t_start=0, t_stop=600)
            np.savez(save_path + '\\' + x + '_spike_data', np.array(curr_spike_data, dtype = np.bool))
        else:
            print("_spike_data file already present")

if __name__ == '__main__':
    data_path = "Z:\\pyDentate\\pyDentateData\\pattern_separation_data_local_input_revised\\seed10003\\scale1000\\net_tunedrev\\"
    for root, dirs, files in os.walk(data_path):
        print(root)
        for name in files:
            if str(name).endswith('.pydd'):
                #try:
                    print("Getting binary_traces")
                    get_binary_traces(root)
                #except:
                    print("Error encountered in get_binary_traces()")
                    break
