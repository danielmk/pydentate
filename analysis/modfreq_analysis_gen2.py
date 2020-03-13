#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 10:34:46 2020

@author: daniel
"""

import numpy as np
import os
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pdb
from scipy import signal

data_path = "/home/daniel/repos/output_dynamics_gen2_long/"

def make_filename(params, run):
    """Generate a filename based on params and run.
    params: (net, seed, freq, weight)"""
    nets = ['net_ppdynamicstuned',
            'net_ppdynamicsoff']

    return(f"time-stamps_{nets[int(params[0])]}.TunedNetwork_{params[2]}_{str(run).zfill(3)}_00.00025_00.00050_{params[1]}_0100_.npz")
#"time-stamps_net_ppdynamicsoff.TunedNetwork_010002_004_00.00025_00.00050_0010_0100_"
files = [x for x in os.listdir(data_path) if ".npz" in x]

params = []
for x in files:
    x_split = x.split('_')
    x_weight = '.'.join(x_split[7].split('.')[0:2])
    if len(x_split) != 10:
        print("ALARM")
        break
    if 'ppdynamicsoff' in x_split[2]:
        net=1
    elif 'ppdynamicstuned' in x_split[2]:
        net=0
    # params: (net, freq, run)
    params.append((net,
                   x_split[7],
                   x_split[3]))

params_array = np.array(params)
params_unique = np.unique(params_array, axis=0)
result_dr = []
result_tdr = []
result_dr_highly_similar = []
result_gc_n_active = []
result_gc_n_spikes = []
result_vec_pp = []
result_vec_gc = []
result_vec_mc = []
result_vec_bc = []
result_vec_hc = []

# Loop through the unique params to calculate results for each run
runs=range(25)
for param in params_unique:
    vec_in_list = []
    vec_out_list = []
    vec_tin_list = []
    vec_tout_list = []
    vec_pp_list = []
    vec_gc_list = []
    vec_mc_list = []
    vec_bc_list = []
    vec_hc_list = []
    for run in runs:
        curr_filename = make_filename(param, run)
        curr_file = np.load(data_path+curr_filename)
        vec_pp_list.append(curr_file['pp_ts'])
        vec_gc_list.append(curr_file['gc_ts'])
        vec_mc_list.append(curr_file['mc_ts'])
        vec_bc_list.append(curr_file['bc_ts'])
        vec_hc_list.append(curr_file['hc_ts'])
        vec_in_list.append(curr_file['pp_ts'][:, 2000:].sum(axis=1))
        vec_out_list.append(curr_file['gc_ts'][:, 2000:].sum(axis=1))
        vec_tin_list.append(curr_file['pp_ts'][:, 2000:].sum(axis=0))
        vec_tout_list.append(curr_file['gc_ts'][:, 2000:].sum(axis=0))
    vec_in = np.stack(vec_in_list, axis=1).transpose()
    vec_out = np.stack(vec_out_list, axis=1).transpose()
    #vec_in_corrs = [pearsonr(x,y)[0] for x in vec_in for y in vec_in]
    vec_in_corrs = [pearsonr(x,y)[0] for x in vec_in for y in vec_in]
    vec_in_corrs = np.array(vec_in_corrs).reshape((25,25))
    vec_in_corrs_triu = np.triu(vec_in_corrs, k=1)
    vec_in_corrs_triu_highly_similar = vec_in_corrs_triu > 0.9
    vec_out_corrs = [pearsonr(x,y)[0] for x in vec_out for y in vec_out]
    vec_out_corrs = np.array(vec_out_corrs).reshape((25,25))
    vec_out_corrs_triu = np.triu(vec_out_corrs, k=1)
    n = 300  # number of elements in the upper triangular 25x25 matrix
    n_highly_similar = vec_in_corrs_triu_highly_similar.sum()
    dr = vec_in_corrs_triu - vec_out_corrs_triu
    dr = np.nanmean(dr[dr!=0])
    dr_highly_similar = (vec_in_corrs_triu[vec_in_corrs_triu_highly_similar] -
                         vec_out_corrs_triu[vec_in_corrs_triu_highly_similar])
    dr_highly_similar = np.nanmean(dr_highly_similar[dr_highly_similar!=0])
    gc_n_active = (vec_out > 0).sum(axis=1).mean()
    gc_avg_spikes = vec_out[vec_out>0].mean()
    result_dr.append(dr)
    result_dr_highly_similar.append(dr_highly_similar)
    result_gc_n_active.append(gc_n_active)
    result_gc_n_spikes.append(gc_avg_spikes)

    vec_tin = np.stack(vec_tin_list, axis=1).transpose()
    vec_tin = vec_tin.reshape((25,int(vec_tin.shape[1]/200),200)).sum(axis=2)
    vec_tout = np.stack(vec_tout_list, axis=1).transpose()
    vec_tout = vec_tout.reshape((25,int(vec_tout.shape[1]/200),200)).sum(axis=2)
    vec_tin_corrs = [pearsonr(x,y)[0] for x in vec_tin for y in vec_tin]
    vec_tin_corrs = np.array(vec_tin_corrs).reshape((25,25))
    vec_tin_corrs_triu = np.triu(vec_tin_corrs, k=1)
    vec_tout_corrs = [pearsonr(x,y)[0] for x in vec_tout for y in vec_tout]
    vec_tout_corrs = np.array(vec_tout_corrs).reshape((25,25))
    vec_tout_corrs_triu = np.triu(vec_tout_corrs, k=1)
    vec_tout_corrs_triu = np.absolute(vec_tout_corrs_triu)
    tdr = (vec_tin_corrs_triu - vec_tout_corrs_triu)
    tdr = np.nanmean(tdr[tdr!=0])
    result_tdr.append(tdr)
    if np.isnan(dr):
        pdb.set_trace()
    result_vec_pp.append(vec_pp_list)
    result_vec_gc.append(vec_gc_list)
    result_vec_mc.append(vec_mc_list)
    result_vec_bc.append(vec_bc_list)
    result_vec_hc.append(vec_hc_list)

result_dr = np.array(result_dr)
result_tdr = np.array(result_tdr)
result_dr_highly_similar = np.array(result_dr_highly_similar)
result_gc_n_active = np.array(result_gc_n_active)
result_gc_n_spikes = np.array(result_gc_n_spikes)
result_vec_pp = np.array(result_vec_pp)
result_vec_gc = np.array(result_vec_gc)
result_vec_mc = np.array(result_vec_mc)
result_vec_bc = np.array(result_vec_bc)
result_vec_hc = np.array(result_vec_hc)

# Now that each run is calculated, we need to average across runs
params_no_runs = np.array([(x[0],x[1]) for x in params_unique])
params_no_runs_unique = np.unique(params_no_runs, axis=0)
params_augmented = np.zeros((params_no_runs_unique.shape[0], params_no_runs.shape[1]+14))

# Go through the unique parameters without runs and collect all runs
for idx1, p1 in enumerate(params_no_runs_unique):
    curr_idc = []
    for idx2, p2 in enumerate(params_unique):
        if p1[0] == p2[0] and p1[1] == p2[1]:
            curr_idc.append(idx2)
    # Check number of runs.
    if len(curr_idc) != 30:
        raise ValueError("Critical problem with number of runs")

    avg_dr = result_dr[curr_idc].mean()
    sem_dr = result_dr[curr_idc].std() / np.sqrt(len(curr_idc))
    avg_tdr = result_tdr[curr_idc].mean()
    sem_tdr = result_tdr[curr_idc].std() / np.sqrt(len(curr_idc))
    avg_gc_n_active = result_gc_n_active[curr_idc].mean()
    sem_gc_n_active = result_gc_n_active[curr_idc].std() / np.sqrt(len(curr_idc))
    avg_gc_n_spikes = result_gc_n_spikes[curr_idc].mean()
    sem_gc_n_spikes = result_gc_n_spikes[curr_idc].std() / np.sqrt(len(curr_idc))
    avg_gc_perc_active = ((result_gc_n_active[curr_idc] / 2000)*100).mean()
    sem_gc_perc_active = ((result_gc_n_active[curr_idc] / 2000)*100).std() / np.sqrt(len(curr_idc))
    avg_dr_highly_similar = result_dr_highly_similar[curr_idc].mean()
    sem_dr_highly_similar = result_dr_highly_similar[curr_idc].std() / np.sqrt(len(curr_idc))
    params_augmented[idx1, 0] = float(p1[0])
    params_augmented[idx1, 1] = float(p1[1])
    params_augmented[idx1, 3] = avg_dr
    params_augmented[idx1, 4] = sem_dr
    params_augmented[idx1, 5] = avg_gc_n_active
    params_augmented[idx1, 6] = sem_gc_n_active
    params_augmented[idx1, 7] = avg_gc_n_spikes
    params_augmented[idx1, 8] = sem_gc_n_spikes
    params_augmented[idx1, 9] = avg_gc_perc_active
    params_augmented[idx1, 10] = sem_gc_perc_active
    params_augmented[idx1, 11] = avg_dr_highly_similar
    params_augmented[idx1, 12] = sem_dr_highly_similar
    params_augmented[idx1, 13] = avg_tdr
    params_augmented[idx1, 14] = sem_tdr

tuned = np.array([x for x in params_augmented if x[0] == 0])
off = np.array([x for x in params_augmented if x[0] == 1])
#tuned_minus_off = 

# Start plotting stuff
ctuned = 'g'
coff = 'b'
legend = (('tuned_10Hz', 'tuned_30Hz', 'nofb_10Hz', 'nofb_30Hz', 'disinh_10Hz', 'disinh_30Hz'))
capsize=5
capthick=1

"""PLOTTING REPRESENTATIVE EXAMPLES"""
# result_vec_gc[1,0,:,:] looks ok
# result_vec_gc[4,0,:,:] looks ok
# result_vec_gc[6,0,:,:] looks ok
idx_tuned = 6
idx_off = 66
time = np.arange(0, 1100, 0.1)
#time2d_pp = np.tile(time1d, 400).reshape((400,11000))
#time2d_gc = np.tile(time1d, 2000).reshape((2000,11000))
pp_tuned_ts = [time[x] for x in result_vec_pp[idx_tuned,0,:,:]]
gc_tuned_ts = [time[x] for x in result_vec_gc[idx_tuned,0,:,:]]
pp_off_ts = [time[x] for x in result_vec_pp[idx_off,0,:,:]]
gc_off_ts = [time[x] for x in result_vec_gc[idx_off,0,:,:]]
fig1, axes = plt.subplots(2,2)
axes[0][0].eventplot(pp_tuned_ts[0:24], color='k')
axes[1][0].eventplot(gc_tuned_ts, color='k', linelengths=8)
axes[0][1].eventplot(pp_off_ts[0:24], color='k')
axes[1][1].eventplot(gc_off_ts, color='k', linelengths=8)
for x in axes:
    for y in x:
        y.set_xlim((-100,1100))
# Plot PP strength against dR
fig1, axes = plt.subplots(4)
"""
vecs=tuned
axes[0].errorbar(x=vecs[:,1], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=ctuned)
axes[1].errorbar(x=vecs[:,1], y=vecs[:,3], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ctuned)
axes[2].errorbar(x=vecs[:,1], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ctuned)
axes[3].errorbar(x=vecs[:,1], y=vecs[:,13], yerr= vecs[:,14], capsize=capsize, capthick=capthick, color=ctuned)
"""
vecs=tuned
axes[0].errorbar(x=vecs[:,1], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=ctuned)
axes[1].errorbar(x=vecs[:,1], y=vecs[:,3], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ctuned)
axes[2].errorbar(x=vecs[:,1], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ctuned)
axes[3].errorbar(x=vecs[:,1], y=vecs[:,13], yerr= vecs[:,14], capsize=capsize, capthick=capthick, color=ctuned)
axes[0].set_xlabel('PP weight')
axes[0].set_ylabel('% Active Cells')
axes[1].set_ylabel('dR')
axes[2].set_ylabel('dR highly similar')
axes[3].set_ylabel('dtR highly similar')

for x in axes:
    x.set_xticks(np.arange(5,101,5))
    
    
# Start plotting stuff
ctuned = 'g'
coff = 'b'
legend = (('tuned_10Hz', 'tuned_30Hz', 'nofb_10Hz', 'nofb_30Hz', 'disinh_10Hz', 'disinh_30Hz'))
capsize=5
capthick=1
