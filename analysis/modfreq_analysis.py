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

data_path = "/home/daniel/repos/output_dynamics_gen2_long/"

def make_filename(params, run):
    """Generate a filename based on params and run.
    params: (net, seed, freq, weight)"""
    nets = ['net_tuneddynamics',
            'net_ppdynamicsoff']

    return(f"time-stamps_{nets[int(params[0])]}.TunedNetwork_{params[2]}_{str(run).zfill(3)}_02000_00060"
           f"_00024_00024_00024_0020_0020_00.00100_00.00100_{params[1]}_0100_00.02500"
           f"_00.02500_00.00120_00.00600.npz")

files = [x for x in os.listdir(data_path) if ".npz" in x]

params = []
for x in files:
    x_split = x.split('_')
    x_weight = '.'.join(x_split[7].split('.')[0:2])
    if len(x_split) != 20:
        print("ALARM")
        break
    if 'ppdynamicsoff' in x_split[2]:
        net=1
    elif 'tuneddynamics' in x_split[2]:
        net=0
    # params: (net, freq, run)
    params.append((net,
                   x_split[14],
                   x_split[3]))

params_array = np.array(params)
params_unique = np.unique(params_array, axis=0)
result_dr = []
result_tdr = []
result_dr_highly_similar = []
result_gc_n_active = []
result_gc_n_spikes = []

# Loop through the unique params to calculate results for each run
runs=range(25)
for param in params_unique:
    vec_in_list = []
    vec_out_list = []
    vec_tin_list = []
    vec_tout_list = []
    for run in runs:
        curr_filename = make_filename(param, run)
        curr_file = np.load(data_path+curr_filename)
        vec_in_list.append(curr_file['pp_ts'][:,1000:].sum(axis=1))
        vec_out_list.append(curr_file['gc_ts'][:,1000:].sum(axis=1))
        vec_tin_list.append(curr_file['pp_ts'][:,1000:].sum(axis=0))
        vec_tout_list.append(curr_file['gc_ts'][:,1000:].sum(axis=0))
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
    dr = (vec_in_corrs_triu - vec_out_corrs_triu).sum()/n
    dr_highly_similar = (vec_in_corrs_triu[vec_in_corrs_triu_highly_similar] -
                         vec_out_corrs_triu[vec_in_corrs_triu_highly_similar]).sum()/n_highly_similar
    gc_n_active = (vec_out > 0).sum(axis=1).mean()
    gc_avg_spikes = vec_out[vec_out>0].mean()
    result_dr.append(dr)
    result_dr_highly_similar.append(dr_highly_similar)
    result_gc_n_active.append(gc_n_active)
    result_gc_n_spikes.append(gc_avg_spikes)

    vec_tin = np.stack(vec_tin_list, axis=1).transpose()
    vec_tin = vec_tin.reshape((25,25,200)).sum(axis=2)
    vec_tout = np.stack(vec_tout_list, axis=1).transpose()
    vec_tout = vec_tout.reshape((25,25,200)).sum(axis=2)
    vec_tin_corrs = [pearsonr(x,y)[0] for x in vec_tin for y in vec_tin]
    vec_tin_corrs = np.array(vec_tin_corrs).reshape((25,25))
    vec_tin_corrs_triu = np.triu(vec_tin_corrs, k=1)
    vec_tout_corrs = [pearsonr(x,y)[0] for x in vec_tout for y in vec_tout]
    vec_tout_corrs = np.array(vec_tout_corrs).reshape((25,25))
    vec_tout_corrs_triu = np.triu(vec_tout_corrs, k=1)
    vec_tout_corrs_triu = np.absolute(vec_tout_corrs_triu)
    tdr = (vec_tin_corrs_triu - vec_tout_corrs_triu).sum()/n
    result_tdr.append(tdr)
    # pdb.set_trace()

result_dr = np.array(result_dr)
result_tdr = np.array(result_tdr)
result_dr_highly_similar = np.array(result_dr_highly_similar)
result_gc_n_active = np.array(result_gc_n_active)
result_gc_n_spikes = np.array(result_gc_n_spikes)

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
    if len(curr_idc) != 10:
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

# Plot PP strength against dR
fig1, axes = plt.subplots(4)
vecs=tuned
axes[0].errorbar(x=vecs[:,1], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=ctuned)
axes[1].errorbar(x=vecs[:,1], y=vecs[:,3], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ctuned)
axes[2].errorbar(x=vecs[:,1], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ctuned)
axes[3].errorbar(x=vecs[:,1], y=vecs[:,13], yerr= vecs[:,14], capsize=capsize, capthick=capthick, color=ctuned)

vecs=off
axes[0].errorbar(x=vecs[:,1], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=coff)
axes[1].errorbar(x=vecs[:,1], y=vecs[:,3], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=coff)
axes[2].errorbar(x=vecs[:,1], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=coff)
axes[3].errorbar(x=vecs[:,1], y=vecs[:,13], yerr= vecs[:,14], capsize=capsize, capthick=capthick, color=coff)
axes[0].set_xlabel('PP weight')
axes[0].set_ylabel('% Active Cells')
axes[1].set_ylabel('dR')
axes[2].set_ylabel('dR highly similar')
axes[3].set_ylabel('dtR highly similar')

for x in axes:
    x.set_xticks(np.arange(5,101,5))

"""
print("Mean +- SEM\n"
      "-----------\n"
      "pp_weight: " + str(tuned_10Hz[:,2]) + '\n'
      "mean_tuned_10Hz: " + str(tuned_10Hz[:,3]) + '\n'
      " SEM_tuned_10Hz: " + str(tuned_10Hz[:,4]) +
      "mean_tuned_30Hz: " + str(tuned_30Hz[:,3]) + '\n'
      " SEM_tuned_30Hz: " + str(tuned_30Hz[:,4]) +
      "mean_nofb_10Hz: " + str(nofb_10Hz[:,3]) + '\n'
      " SEM_nofb_10Hz: " + str(nofb_10Hz[:,4]) + 
      "mean_nofb_30Hz: " + str(nofb_30Hz[:,3]) + '\n'
      " SEM_nofb_30Hz: " + str(nofb_30Hz[:,4]) +
      "mean_disinh_10Hz: " + str(disinh_10Hz[:,3]) + '\n'+
      " SEM_disinh_10Hz: " + str(disinh_10Hz[:,4]) +
      "mean_disinh_30Hz: " + str(disinh_30Hz[:,3]) + '\n'
      " SEM_disinh_30Hz: " + str(disinh_30Hz[:,4]))

vecs=tuned_10Hz
axes[1].errorbar(x=vecs[:,2], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes[1].errorbar(x=vecs[:,2], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes[1].errorbar(x=vecs[:,2], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=cn10Hz)
vecs=nofb_30Hz
axes[1].errorbar(x=vecs[:,2], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=cn30Hz)
vecs=disinh_10Hz
axes[1].errorbar(x=vecs[:,2], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=cd10Hz)
vecs=disinh_30Hz
axes[1].errorbar(x=vecs[:,2], y=vecs[:,9], yerr= vecs[:,10], capsize=capsize, capthick=capthick, color=cd30Hz)
axes[1].set_xlabel('PP weight')
axes[1].set_ylabel('% Active GCs')
print("Mean +- SEM\n"
      "-----------\n"
      "pp_weight: " + str(tuned_10Hz[:,2]) + '\n'
      "mean_tuned_10Hz: " + str(tuned_10Hz[:,9]) + '\n'
      " SEM_tuned_10Hz: " + str(tuned_10Hz[:,10]) +
      "mean_tuned_30Hz: " + str(tuned_30Hz[:,9]) + '\n'
      " SEM_tuned_30Hz: " + str(tuned_30Hz[:,10]) +
      "mean_nofb_10Hz: " + str(nofb_10Hz[:,9]) + '\n'
      " SEM_nofb_10Hz: " + str(nofb_10Hz[:,10]) + 
      "mean_nofb_30Hz: " + str(nofb_30Hz[:,9]) + '\n'
      " SEM_nofb_30Hz: " + str(nofb_30Hz[:,10]) +
      "mean_disinh_10Hz: " + str(disinh_10Hz[:,9]) + '\n'
      " SEM_disinh_10Hz: " + str(disinh_10Hz[:,10]) +
      "mean_disinh_30Hz: " + str(disinh_30Hz[:,9]) + '\n'
      " SEM_disinh_30Hz: " + str(disinh_30Hz[:,10]))

vecs=tuned_10Hz
axes[2].errorbar(x=vecs[:,2], y=vecs[:,7], yerr= vecs[:,8], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes[2].errorbar(x=vecs[:,2], y=vecs[:,7], yerr= vecs[:,8], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes[2].errorbar(x=vecs[:,2], y=vecs[:,7], yerr= vecs[:,8], capsize=capsize, capthick=capthick, color=cn10Hz)
vecs=nofb_30Hz
axes[2].errorbar(x=vecs[:,2], y=vecs[:,7], yerr= vecs[:,8], capsize=capsize, capthick=capthick, color=cn30Hz)
vecs=disinh_10Hz
axes[2].errorbar(x=vecs[:,2], y=vecs[:,7], yerr= vecs[:,8], capsize=capsize, capthick=capthick, color=cd10Hz)
vecs=disinh_30Hz
axes[2].errorbar(x=vecs[:,2], y=vecs[:,7], yerr= vecs[:,8], capsize=capsize, capthick=capthick, color=cd30Hz)
axes[2].set_xlabel('PP weight')
axes[2].set_ylabel('APs/active')
axes[2].legend(legend)

fig2, axes2 = plt.subplots(1,1)
vecs=tuned_10Hz
axes2.errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes2.errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes2.errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cn10Hz)
vecs=nofb_30Hz
axes2.errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cn30Hz)
vecs=disinh_10Hz
axes2.errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cd10Hz)
vecs=disinh_30Hz
axes2.errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cd30Hz)
axes2.legend((legend))
axes2.set_xlabel('% Active GCs')
axes2.set_ylabel('dR')

fig3, axes2 = plt.subplots(3,1)
vecs=tuned_10Hz
axes2[0].errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes2[0].errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes2[1].errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=nofb_30Hz
axes2[1].errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=disinh_10Hz
axes2[2].errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=disinh_30Hz
axes2[2].errorbar(x=vecs[:,9], y=vecs[:,3], xerr=vecs[:,10], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
axes2[2].legend((legend))
axes2[2].set_xlabel('% Active GCs')
axes2[2].set_ylabel('dR')
for x in axes2:
    x.set_xlim((0,40))
    x.set_ylim((0,0.35))

fig4, axes3 = plt.subplots(1,1)
vecs=tuned_10Hz
axes3.errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes3.errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes3.errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cn10Hz)
vecs=nofb_30Hz
axes3.errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cn30Hz)
vecs=disinh_10Hz
axes3.errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cd10Hz)
vecs=disinh_30Hz
axes3.errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=cd30Hz)
axes3.legend((legend))
axes3.set_xlabel('APs/active')
axes3.set_ylabel('dR')

fig5, axes3 = plt.subplots(3,1)
vecs=tuned_10Hz
axes3[0].errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes3[0].errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes3[1].errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=nofb_30Hz
axes3[1].errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=disinh_10Hz
axes3[2].errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=disinh_30Hz
axes3[2].errorbar(x=vecs[:,7], y=vecs[:,3], xerr=vecs[:,8], yerr= vecs[:,4], capsize=capsize, capthick=capthick, color=ct30Hz)
axes3[2].legend((legend))
axes3[2].set_xlabel('APs/active')
axes3[2].set_ylabel('dR')
"""
"""FIGURE S7 STYLE COMPARISON WITHIN INPUT STRENGTH
params_augmented[idx1, 0] = float(p1[0])
params_augmented[idx1, 1] = float(p1[1])
params_augmented[idx1, 2] = float(p1[2])
params_augmented[idx1, 3] = avg_dr
params_augmented[idx1, 4] = sem_dr
params_augmented[idx1, 5] = avg_gc_n_active
params_augmented[idx1, 6] = sem_gc_n_active
params_augmented[idx1, 7] = avg_gc_n_spikes
params_augmented[idx1, 8] = sem_gc_n_spikes
params_augmented[idx1, 9] = avg_gc_perc_active
params_augmented[idx1, 10] = sem_gc_perc_active
"""
"""
fig6, axes4 = plt.subplots(8, 2)
markersize=3

inputs =np.arange(0.0006, 0.002, 0.0002)
for idx, inp in enumerate(inputs):
    x = np.array([tuned_10Hz[idx, 9], tuned_30Hz[idx, 9], nofb_10Hz[idx, 9], nofb_30Hz[idx, 9], disinh_10Hz[idx,9], disinh_30Hz[idx,9]])
    #xerr = np.array([tuned_10Hz[idx, 10], tuned_30Hz[idx, 10], nofb_10Hz[idx, 10], nofb_30Hz[idx, 10], disinh_10Hz[idx,10], disinh_30Hz[idx,10]])
    y = np.array([tuned_10Hz[idx, 3], tuned_30Hz[idx, 3], nofb_10Hz[idx, 3], nofb_30Hz[idx, 3], disinh_10Hz[idx, 3], disinh_30Hz[idx,3]])
    yerr = np.array([tuned_10Hz[idx, 4], tuned_30Hz[idx, 4], nofb_10Hz[idx, 4], nofb_30Hz[idx, 4], disinh_10Hz[idx, 4], disinh_30Hz[idx,4]])
    axes4[idx,0].errorbar(x[0::2], y[0::2], yerr=yerr[0::2], marker='o',linestyle='none', color='blue', markersize=markersize)
    axes4[idx,0].errorbar(x[1::2], y[1::2], yerr=yerr[1::2], marker='o',linestyle='none', color='green', markersize=markersize)
    linreg_10Hz = linregress(x[0::2], y[0::2])
    linreg_30Hz = linregress(x[1::2], y[1::2])
    axes4[idx,0].plot([0,100], [linreg_10Hz[1], linreg_10Hz[1] + linreg_10Hz[0]*100], color='blue', linestyle='-')
    axes4[idx,0].plot([0,100], [linreg_30Hz[1], linreg_30Hz[1] + linreg_30Hz[0]*100], color='green', linestyle='-')
    axes4[idx,0].set_xlim((x.min()-0.05, x.max()+0.05))
    #axes4[idx].set_ylim((y.min()-0.01, y.max()+0.01))
    axes4[idx,0].set_ylim((0, 0.4))
    print(str(inp) + " perc active, 10Hz slope + intercept:\n" + str(linreg_10Hz[0]) +'+' + str(linreg_10Hz[1]) + '\n')
    print(str(inp) + " perc active, 30Hz slope + intercept:\n" + str(linreg_30Hz[0]) +'+' + str(linreg_30Hz[1])+ '\n')
    
for idx, inp in enumerate(inputs):
    x = np.array([tuned_10Hz[idx, 7], tuned_30Hz[idx, 7], nofb_10Hz[idx, 7], nofb_30Hz[idx, 7], disinh_10Hz[idx,7], disinh_30Hz[idx,7]])
    #xerr = np.array([tuned_10Hz[idx, 8], tuned_30Hz[idx, 8], nofb_10Hz[idx, 8], nofb_30Hz[idx, 8], disinh_10Hz[idx,8], disinh_30Hz[idx,8]])
    y = np.array([tuned_10Hz[idx, 3], tuned_30Hz[idx, 3], nofb_10Hz[idx, 3], nofb_30Hz[idx, 3], disinh_10Hz[idx, 3], disinh_30Hz[idx,3]])
    yerr = np.array([tuned_10Hz[idx, 4], tuned_30Hz[idx, 4], nofb_10Hz[idx, 4], nofb_30Hz[idx, 4], disinh_10Hz[idx, 4], disinh_30Hz[idx,4]])
    axes4[idx,1].errorbar(x[0::2], y[0::2], yerr=yerr[0::2], marker='o',linestyle='none', color='blue', markersize=markersize)
    axes4[idx,1].errorbar(x[1::2], y[1::2], yerr=yerr[1::2], marker='o',linestyle='none', color='green', markersize=markersize)
    linreg_10Hz = linregress(x[0::2], y[0::2])
    linreg_30Hz = linregress(x[1::2], y[1::2])
    axes4[idx,1].plot([0,100], [linreg_10Hz[1], linreg_10Hz[1] + linreg_10Hz[0]*100], color='blue', linestyle='-')
    axes4[idx,1].plot([0,100], [linreg_30Hz[1], linreg_30Hz[1] + linreg_30Hz[0]*100], color='green', linestyle='-')
    axes4[idx,1].set_xlim((x.min()-0.05, x.max()+0.05))
    #axes4[idx].set_ylim((y.min()-0.01, y.max()+0.01))
    axes4[idx,1].set_ylim((0, 0.4))
    print(str(inp) + " num APs, 10Hz slope + intercept:\n" + str(linreg_10Hz[0]) +' + ' + str(linreg_10Hz[1])+ '\n')
    print(str(inp) + " num APs, 30Hz slope + intercept:\n" + str(linreg_30Hz[0]) +' + ' + str(linreg_30Hz[1])+ '\n')

# Plot PP strength against dR
fig6, axes6 = plt.subplots(1)
vecs=tuned_10Hz
axes6.errorbar(x=vecs[:,2], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes6.errorbar(x=vecs[:,2], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes6.errorbar(x=vecs[:,2], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=cn10Hz)
vecs=nofb_30Hz
axes6.errorbar(x=vecs[:,2], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=cn30Hz)
vecs=disinh_10Hz
axes6.errorbar(x=vecs[:,2], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=cd10Hz)
vecs=disinh_30Hz
axes6.errorbar(x=vecs[:,2], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=cd30Hz)
axes6.set_xlabel('PP weight')
axes6.set_ylabel('dR')
axes6.legend(legend)
print("Mean +- SEM\n"
      "-----------\n"
      "pp_weight: " + str(tuned_10Hz[:,2]) + '\n'
      "mean_tuned_10Hz: " + str(tuned_10Hz[:,3]) + '\n'
      " SEM_tuned_10Hz: " + str(tuned_10Hz[:,4]) +
      "mean_tuned_30Hz: " + str(tuned_30Hz[:,3]) + '\n'
      " SEM_tuned_30Hz: " + str(tuned_30Hz[:,4]) +
      "mean_nofb_10Hz: " + str(nofb_10Hz[:,3]) + '\n'
      " SEM_nofb_10Hz: " + str(nofb_10Hz[:,4]) + 
      "mean_nofb_30Hz: " + str(nofb_30Hz[:,3]) + '\n'
      " SEM_nofb_30Hz: " + str(nofb_30Hz[:,4]) +
      "mean_disinh_10Hz: " + str(disinh_10Hz[:,3]) + '\n'+
      " SEM_disinh_10Hz: " + str(disinh_10Hz[:,4]) +
      "mean_disinh_30Hz: " + str(disinh_30Hz[:,3]) + '\n'
      " SEM_disinh_30Hz: " + str(disinh_30Hz[:,4]))

fig7, axes7 = plt.subplots(3,1)
vecs=tuned_10Hz
axes7[0].errorbar(x=vecs[:,9], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=tuned_30Hz
axes7[0].errorbar(x=vecs[:,9], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=nofb_10Hz
axes7[1].errorbar(x=vecs[:,9], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=nofb_30Hz
axes7[1].errorbar(x=vecs[:,9], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct30Hz)
vecs=disinh_10Hz
axes7[2].errorbar(x=vecs[:,9], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct10Hz)
vecs=disinh_30Hz
axes7[2].errorbar(x=vecs[:,9], y=vecs[:,11], yerr= vecs[:,12], capsize=capsize, capthick=capthick, color=ct30Hz)
axes7[2].legend((legend))
axes7[2].set_xlabel('% Active GCs')
axes7[2].set_ylabel('dR')
for x in axes7:
    x.set_xlim((0,40))
    x.set_ylim((0,0.2))

print("Mean +- SEM\n"
      "-----------\n"
      "pp_weight: " + str(tuned_10Hz[:,2]) + '\n'
      "mean_tuned_10Hz: " + str(tuned_10Hz[:,11]) + '\n'
      " SEM_tuned_10Hz: " + str(tuned_10Hz[:,12]) +
      "mean_tuned_30Hz: " + str(tuned_30Hz[:,11]) + '\n'
      " SEM_tuned_30Hz: " + str(tuned_30Hz[:,12]) +
      "mean_nofb_10Hz: " + str(nofb_10Hz[:,11]) + '\n'
      " SEM_nofb_10Hz: " + str(nofb_10Hz[:,12]) + 
      "mean_nofb_30Hz: " + str(nofb_30Hz[:,11]) + '\n'
      " SEM_nofb_30Hz: " + str(nofb_30Hz[:,12]) +
      "mean_disinh_10Hz: " + str(disinh_10Hz[:,11]) + '\n'+
      " SEM_disinh_10Hz: " + str(disinh_10Hz[:,12]) +
      "mean_disinh_30Hz: " + str(disinh_30Hz[:,11]) + '\n'
      " SEM_disinh_30Hz: " + str(disinh_30Hz[:,12]))
"""