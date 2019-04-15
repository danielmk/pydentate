# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 10:25:41 2018

@author: daniel
"""
import shelve
import numpy as np
import matplotlib.pyplot as plt
import os

#Home PC
#directory = "C:\\Users\\daniel\\repos\\pyDentate\paradigm_pattern-separation_saves_2018-03-11\\"
#Office PC
#directory = "Y:\\DanielM\\023_Dentate Gyrus Model\\paradigm_spatial-inhibition\\"
#Dropbox
def get_spike_stats(data_path):
    data_files = [f for f in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, f)) and '.pydd' in f and not '.npz' in f]
    data_files.sort()
    
    perc_active_gcs_list = []
    perc_active_mcs_list = []
    perc_active_bcs_list = []
    perc_active_hcs_list = []
    
    n_aps_avg_gcs_list = []
    n_aps_avg_mcs_list = []
    n_aps_avg_bcs_list = []
    n_aps_avg_hcs_list = []
    
    n_aps_std_gcs_list = []
    n_aps_std_mcs_list = []
    n_aps_std_bcs_list = []
    n_aps_std_hcs_list = []
    
    # Get to BasketCell Connection
    for x in data_files:
        curr_data = shelve.open(data_path + x)
        active_gcs = np.array(np.argwhere(np.array(curr_data[curr_data.keys()[0]]['populations'][0]['ap_number']) > 0),dtype = int).flatten()
        active_mcs = np.array(np.argwhere(np.array(curr_data[curr_data.keys()[0]]['populations'][1]['ap_number']) > 0),dtype = int).flatten()
        active_bcs = np.array(np.argwhere(np.array(curr_data[curr_data.keys()[0]]['populations'][2]['ap_number']) > 0),dtype = int).flatten()
        active_hcs = np.array(np.argwhere(np.array(curr_data[curr_data.keys()[0]]['populations'][3]['ap_number']) > 0),dtype = int).flatten()
    
        n_aps_avg_gcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][0]['ap_number'])[active_gcs].mean())
        n_aps_std_gcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][0]['ap_number'])[active_gcs].std())
    
        n_aps_avg_mcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][1]['ap_number'])[active_mcs].mean())
        n_aps_std_mcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][1]['ap_number'])[active_mcs].std())
    
        n_aps_avg_bcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][2]['ap_number'])[active_bcs].mean())
        n_aps_std_bcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][2]['ap_number'])[active_bcs].std())
        
        n_aps_avg_hcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][3]['ap_number'])[active_hcs].mean())
        n_aps_std_hcs_list.append(np.array(curr_data[curr_data.keys()[0]]['populations'][3]['ap_number'])[active_hcs].std())
    
        perc_active_gcs_list.append((len(active_gcs) / 2000.0)*100)
        perc_active_mcs_list.append((len(active_mcs) / 60.0)*100)
        perc_active_bcs_list.append((len(active_bcs) / 24.0)*100)
        perc_active_hcs_list.append((len(active_hcs) / 24.0)*100)
    
    perc_active_cells = np.array((perc_active_gcs_list, perc_active_mcs_list, perc_active_bcs_list, perc_active_hcs_list))
    n_aps_avg = np.array((n_aps_avg_gcs_list, n_aps_avg_mcs_list, n_aps_avg_bcs_list, n_aps_avg_hcs_list))
    n_aps_std = np.array((n_aps_std_gcs_list, n_aps_std_mcs_list, n_aps_std_bcs_list, n_aps_std_hcs_list))
    
    np.savetxt(data_path + "1_perc_active_cells.txt", perc_active_cells, delimiter = '\t')
    np.savetxt(data_path + "1_n_aps_avg.txt", n_aps_avg, delimiter = '\t')
    np.savetxt(data_path + "1_n_aps_std.txt", n_aps_std, delimiter = '\t')

if __name__ == '__main__':
    data_path = "C:\\Users\\Daniel\\pyDentateData\\pattern_separation_data_local_input_revised\\seed10000\\scale500\\net_nofeedbackrev\\"
    get_spike_stats(data_path)

