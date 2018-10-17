# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
# from scipy.stats import binned_statistic
from scipy.stats import wilcoxon
import pandas as pd
# from scipy.optimize import curve_fit

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# what to analyse
basePath = Path(r"R:\pyDentate\pyDentateData")
model_input = r"pattern_separation_data_local_input_revised"
scale ="scale1000"
filename = r"1_leutgeb-measure_matrix_len-bin_6000.txt"

dictionary = {}
seeds = ['seed10000', 'seed10001', 'seed10002', 'seed10003',
               'seed10004', 'seed10005', 'seed10006']

for seed in seeds:
    dictionary[seed] = {}
    case = r'input_patterns'
     # inputs
    file_to_open = basePath / model_input / seed / case / filename
    try:
        data = np.loadtxt(file_to_open)
        linear_data = np.array([])
        for x in range(25):
            linear_data = np.concatenate((linear_data,
                                          data[x, x:data.shape[0]]))
        dictionary[seed]['inputs'] = linear_data
    except:
        print("no input for " + seed)

    for case in [r'net_nofeedbackrev', r'net_globalrev',
                 r'net_tunedrev', r'net_nonfacilitatingrev',
                 r'net_disinhibitedrev', r'net_reshuffledrev']:
 
        file_to_open = basePath / model_input / seed / scale / case / filename
        try:
            data = np.loadtxt(file_to_open)
            linear_data = np.array([])
            for x in range(25):
                linear_data = np.concatenate((linear_data,
                                              data[x, x:data.shape[0]]))
            dictionary[seed][case] = linear_data
        except:
            print("no " + case + " in " + seed)
            dictionary[seed][case] = dictionary[seed]['inputs']

def show_scatter():
    f1, ax_array = plt.subplots(7, 5, figsize=(11, 16))
    f1.subplots_adjust(hspace=0.1, wspace=0.1)
    
    results = np.zeros([len(seeds), 5])
    for i, seed in enumerate(seeds):
        inputs = dictionary[seed]['inputs']
        tuned = dictionary[seed][r'net_tunedrev']
        global_i = dictionary[seed][r'net_globalrev']
        disinh = dictionary[seed][r'net_disinhibitedrev']
        noFB = dictionary[seed][r'net_nofeedbackrev']
        nonfac = dictionary[seed][r'net_nonfacilitatingrev']
        
        ax = ax_array[i, 0]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, tuned)
        ax.scatter(inputs, tuned, alpha=0.2, c='b', marker='.')
        ax.set_ylabel('output R \n' + str(seed))
        # ax.set_xlabel('input R')
        ax.plot([0, 1], c='k', linestyle='-')
        ax.plot(bins, means, c='k', linestyle='--')
        if(i==0):
            ax.set_title('tuned FB inh.', fontsize=10)
            ax.tick_params(labelbottom='off')
        elif(i==len(seeds)-1):
            ax.set_xlabel('input R')
        else:
            ax.tick_params(labelbottom='off')
        results[i, 0] = np.mean(bins[1:-1] - means[1:-1]) 
    
        ax = ax_array[i, 1]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, noFB)
        ax.scatter(inputs, noFB, alpha=0.2, c='b', marker='.')
        ax.tick_params(labelleft='off')
        # ax.set_xlabel('input R')
        ax.plot([0, 1], c='k', linestyle='-')
        ax.plot(bins, means, c='k', linestyle='--')
        if(i==0):
            ax.set_title('no FB inh.', fontsize=10)
            ax.tick_params(labelbottom='off')
        elif(i==len(seeds)-1):
            ax.set_xlabel('input R')
        else:
            ax.tick_params(labelbottom='off')
        results[i, 1] = np.mean(bins[1:-1] - means[1:-1]) 
    
        ax = ax_array[i, 2]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, disinh)
        ax.scatter(inputs, disinh, alpha=0.2, c='b', marker='.')
        ax.tick_params(labelleft='off')
        # ax.set_xlabel('input R')
        ax.plot([0, 1], c='k', linestyle='-')
        ax.plot(bins, means, c='k', linestyle='--')
        if(i==0):
            ax.set_title('no inh.', fontsize=10)
            ax.tick_params(labelbottom='off')
        elif(i==len(seeds)-1):
            ax.set_xlabel('input R')
        else:
            ax.tick_params(labelbottom='off')
        results[i, 2] = np.mean(bins[1:-1] - means[1:-1]) 
    
        ax = ax_array[i, 3]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, global_i)
        ax.scatter(inputs, global_i, alpha=0.2, c='b', marker='.')
        ax.tick_params(labelleft='off')
        # ax4.set_xlabel('input R')
        ax.plot([0, 1], c='k', linestyle='-')
        ax.plot(bins, means, c='k', linestyle='--')
        if(i==0):
            ax.set_title('global', fontsize=10)
            ax.tick_params(labelbottom='off')
        elif(i==len(seeds)-1):
            ax.set_xlabel('input R')
        else:
            ax.tick_params(labelbottom='off')
        results[i, 3] = np.mean(bins[1:-1] - means[1:-1]) 
        
        ax = ax_array[i, 4]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, nonfac)
        ax.scatter(inputs, nonfac, alpha=0.2, c='b', marker='.')
        ax.tick_params(labelleft='off')
        # ax5.set_xlabel('input R')
        ax.plot([0, 1], c='k', linestyle='-')
        ax.plot(bins, means, c='k', linestyle='--')
        if(i==0):
            ax.set_title('nonfacilitating', fontsize=10)
            ax.tick_params(labelbottom='off')
        elif(i==len(seeds)-1):
            ax.set_xlabel('input R')
        else:
            ax.tick_params(labelbottom='off')
        results[i, 4] = np.mean(bins[1:-1] - means[1:-1]) 
 
    savename = 'I-O-corr_'+ model_input + scale
    savefigname = savename + '.pdf'
    savetxtname = savename + '.txt'
    f1.savefig(savefigname)
    np.savetxt(savetxtname, results)
    
    f2, ax_array = plt.subplots(1, 5, figsize=(11, 2))
    f2.subplots_adjust(hspace=0.1, wspace=0.4)
    
    ax = ax_array[0]
    x_bar = [1, 2, 3]
    bar1 = results[:, 0]
    bar2 = results[:, 1]
    bar3 = results[:, 2]
    means = [np.mean(bar1), np.mean(bar2), np.mean(bar3)]
    SDs = [np.std(bar1), np.std(bar2), np.std(bar3)]
    SEMs = SDs/np.sqrt(len(bar1))
    ax.bar(x_bar, means, yerr=SEMs, color='b', alpha=0.8)
    ax.set_xticklabels(['','tuned','no FB','no inh.'], rotation= 30, fontsize=10)
    ax.set_ylim([0,0.24])
    for seed in range(len(seeds)):
        ax.plot(x_bar, results[seed,:3], c='k', alpha=0.5)

    ax = ax_array[1]
    x_bar = [1, 2, 3]
    noFB = results[:, 1]
    bar1 = results[:, 3] - results[:, 1]
    bar2 = results[:, 0] - results[:, 1]
    bar3 = results[:, 4] - results[:, 1]
    means = [np.mean(bar1), np.mean(bar2), np.mean(bar3)]
    SDs = [np.std(bar1), np.std(bar2), np.std(bar3)]
    SEMs = SDs/np.sqrt(len(bar1))
    ax.bar(x_bar, means, yerr=SEMs, color='b', alpha=0.8)
    ax.set_xticklabels(['','global','tuned','nonfacilitating'], rotation= 30, fontsize=10)
    ax.set_ylim([0,0.13])
    for seed in range(len(seeds)):
        ax.plot(x_bar, results[seed,[3,0,4]] - noFB[seed], c='k', alpha=0.5)

    savename = 'I-O-corr_'+ model_input + scale
    savefigname = savename + '.pdf'
    f2.savefig(savefigname)
    
    return results

def bindata(inputs, linear_data):
    data = linear_data
    bins = np.linspace(0.01, 0.99, 11)
    digitized = np.digitize(inputs, bins) - 1
    binned_data = [data[digitized == i] for i in range(len(bins))]
    # p = [wilcoxon(binned_data[i][:])[1] for i in range(len(bins))]
    means = np.array([data[digitized == i].mean() for i in range(len(bins))])
    SDs = [data[digitized == i].std() for i in range(len(bins))]
    n = [len(binned_data[i][:]) for i in range(len(bins))]
    SEMs = SDs/np.sqrt(n)
    significant = []
    marks = []
    # plt.plot(x, means, color='b')
    # plt.fill_between(x, means - SDs, means + SDs, alpha=0.2, color='b')
        # plt.scatter(significant, marks, color='k', marker='.')
    # significant = bins[np.array(p) < 0.05/10]
    # marks = np.ones(len(significant))
    means = np.append(0, means)
    SDs = np.append(0, SDs)
    SEMs = np.append(0, SEMs)
    bins = bins + 0.05  # bin middles for plotting
    bins[-1] = 1
    bins = np.append(0, bins)

    return significant, marks, means, SDs, SEMs, bins


results = show_scatter()