#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
paradigm_drives gives one .npz file for each run. This script merges
those files into one file for each simulation condition.

@author: daniel
"""

import numpy as np
import os
import matplotlib.pyplot as plt

load_path = "/home/daniel/repos/output/"
save_path = "/home/daniel/repos/pyDentate2Data/"

ts_files = [x for x in os.listdir(load_path) if x[:11] == 'time-stamps']
first_run_files = [x for x in ts_files if x.split('_')[3] == '010000' and x.split('_')[4] == '000']

master_list = []
for x in first_run_files:
    condition_signature = x[58:]
    curr_list = []
    for y in ts_files:
        if y[58:] == condition_signature:
            curr_list.append(y)
    if len(curr_list) == 25:
        curr_pp_list = []
        curr_gc_list = []
        curr_mc_list = []
        curr_bc_list = []
        curr_hc_list = []
        for z in curr_list:
            curr_file = np.load(load_path+z)
            curr_pp_list.append(curr_file['pp_ts'])
            curr_gc_list.append(curr_file['gc_ts'])
            curr_mc_list.append(curr_file['mc_ts'])
            curr_bc_list.append(curr_file['bc_ts'])
            curr_hc_list.append(curr_file['hc_ts'])
            curr_file.close()

        split_fname = x.split('_')
        new_file_name = '_'.join(split_fname[0:3]) + '_25-runs_' +  '_'.join(split_fname[5:])

        np.savez(save_path+new_file_name,
                 pp_ts = np.array(curr_pp_list),
                 gc_ts = np.array(curr_gc_list),
                 mc_ts = np.array(curr_mc_list),
                 bc_ts = np.array(curr_bc_list),
                 hc_ts = np.array(curr_hc_list))
        for z in curr_list:
            os.rename(load_path+z, load_path+'/merged/'+z)