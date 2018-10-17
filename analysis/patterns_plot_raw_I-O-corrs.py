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
# import pandas as pd
# from scipy.optimize import curve_fit

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# what to analyse
basePath = Path(r"R:\pyDentate\pyDentateData")
model_input = r"pattern_separation_data_local_repeated_input_revised"
seed = "seed10000"
scale = "scale1000"
filename = r"1_leutgeb-measure_matrix_len-bin_6000.txt"


case = r"input_patterns_seed_10000"
file_to_open = basePath / model_input / seed / case / filename
data = np.loadtxt(file_to_open)
linear_data= np.array([])
for x in range(25):
    linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
inputs = linear_data

case = r"net_disinhibitedrev"
file_to_open = basePath / model_input / seed / scale /case / filename
try:
    data = np.loadtxt(file_to_open)
    linear_data= np.array([])
    for x in range(25):
        linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
    disinhibited = linear_data
except:
    disinhibited = inputs
    print("no disinhibited data")

case=r"net_nofeedbackrev"
file_to_open = basePath / model_input / seed / scale /case / filename
try:
    data = np.loadtxt(file_to_open)
    linear_data= np.array([])
    for x in range(25):
        linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
    noFeedBack = linear_data
except:
    noFeedBack = inputs
    print("no noFeedBack data")

case = r"net_tunedrev"
file_to_open = basePath / model_input / seed / scale /case / filename
try:
    data = np.loadtxt(file_to_open)
    linear_data= np.array([])
    for x in range(25):
        linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
    tuned = linear_data
except:
    tuned = inputs
    print("no tuned data")

case = r"net_globalrev"
file_to_open = basePath / model_input / seed / scale /case / filename
try:
    data = np.loadtxt(file_to_open)
    linear_data= np.array([])
    for x in range(25):
        linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
    global_i = linear_data
except:
    global_i = inputs
    print("no global data")
    
case = r"net_nonfacilitatingrev"
file_to_open = basePath / model_input / seed / scale /case / filename
try:
    data = np.loadtxt(file_to_open)
    linear_data= np.array([])
    for x in range(25):
        linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
    nonfacilitating = linear_data
except:
    nonfacilitating = inputs
    print("no nonfacilitating data")

def show_scatter():
    f1, ((ax1, ax2, ax3, ax4, ax5), (ax6, ax7, ax8, ax9, ax10)) = plt.subplots(2, 5, figsize=(11, 5))
    f1.subplots_adjust(hspace=0.2, wspace=0.4)

    ax1.scatter(inputs, tuned, alpha=0.2, c='b', marker='.')
    ax1.set_title('tuned FBI', fontsize=10)
    ax1.set_ylabel('output R')
    #  ax1.set_xlabel('input R')
    ax1.plot([0, 1], c='k')

    ax2.scatter(inputs, global_i, alpha=0.2, c='b', marker='.')
    ax2.set_title('global FBI', fontsize=10)
    # ax2.set_xlabel('input R')
    ax2.plot([0, 1], c='k')

    ax3.scatter(inputs, nonfacilitating, alpha=0.2, c='b', marker='.')
    ax3.set_title('nonfac. FBI', fontsize=10)
    # ax3.set_xlabel('input R')
    ax3.plot([0, 1], c='k')

    ax4.scatter(inputs, noFeedBack, alpha=0.2, c='b', marker='.')
    ax4.set_title('no FBI', fontsize=10)
    # ax4.set_xlabel('input R')
    ax4.plot([0, 1], c='k')

    ax5.scatter(inputs, disinhibited, alpha=0.2, c='b', marker='.')
    ax5.set_title('no inh.', fontsize=10)
    # ax5.set_xlabel('input R')
    ax5.plot([0, 1], c='k')

    significant, marks, means, SDs, bins = bindata(inputs, tuned)
    ax6.plot(bins, means, color='b')
    ax6.fill_between(bins, means - SDs, means + SDs, alpha=0.2, color='b')
    ax6.set_ylabel('output R')
    ax6.set_xlabel('input R')
    ax6.plot([0, 1], c='k')

    significant, marks, means, SDs, bins = bindata(inputs, global_i)
    ax7.plot(bins, means, color='b')
    ax7.fill_between(bins, means - SDs, means + SDs, alpha=0.2, color='b')
    ax7.set_ylabel('output R')
    ax7.set_xlabel('input R')
    ax7.plot([0, 1], c='k')
    
    significant, marks, means, SDs, bins = bindata(inputs, nonfacilitating)
    ax8.plot(bins, means, color='b')
    ax8.fill_between(bins, means - SDs, means + SDs, alpha=0.2, color='b')
    ax8.set_ylabel('output R')
    ax8.set_xlabel('input R')
    ax8.plot([0, 1], c='k')
    
    significant, marks, means, SDs, bins = bindata(inputs, noFeedBack)
    ax9.plot(bins, means, color='b')
    ax9.fill_between(bins, means - SDs, means + SDs, alpha=0.2, color='b')
    ax9.set_ylabel('output R')
    ax9.set_xlabel('input R')
    ax9.plot([0, 1], c='k')
    
    significant, marks, means, SDs, bins = bindata(inputs, disinhibited)
    ax10.plot(bins, means, color='b')
    ax10.fill_between(bins, means - SDs, means + SDs, alpha=0.2, color='b')
    ax10.set_ylabel('output R')
    ax10.set_xlabel('input R')
    ax10.plot([0, 1], c='k')
    
    savename = 'I-O-corr_'+model_input + seed + scale + '.pdf'
    f1.savefig(savename)

def bindata(inputs, linear_data):
    data = linear_data
    bins = np.linspace(0.05, 0.95, 10)
    digitized = np.digitize(inputs, bins)
    binned_data = [data[digitized == i] for i in range(len(bins))]
    p = [wilcoxon(binned_data[i][:])[1] for i in range(len(bins))]
    means = np.array([data[digitized == i].mean() for i in range(len(bins))])
    SDs = [data[digitized == i].std() for i in range(len(bins))]
    
    significant = bins[np.array(p)<0.05/10]
    marks = np.ones(len(significant))*0.2
    means = np.append(np.append(0, means), 1)
    SDs = np.append(np.append(0, SDs), 0)
    bins = np.append(np.append(0, bins), 1)
    # plt.plot(bins, means, color='b')
    # plt.fill_between(bins, means - SDs, means + SDs, alpha=0.2, color='b')
    # plt.scatter(significant, marks, color='k', marker='.')
    return significant, marks, means, SDs, bins
    

def show_diff():
    f1, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(11, 2.5))
    f1.subplots_adjust(hspace=0.5,wspace=0.4)
    
    ax1.scatter(inputs,disinhibited-noFeedBack, alpha=0.2, c='b', marker='.')
    ax1.set_ylabel('no feedback inhibition')
    ax1.set_xlabel('input R')
    ax1.set_ylim([-0.02,0.23])

    ax2.scatter(inputs,disinhibited-tuned, alpha=0.2, c='b', marker='.')
    ax2.set_ylabel('tuned feedback inh.')
    ax2.set_xlabel('input R')
    ax2.set_ylim([-0.02,0.23])
    
    ax3.scatter(inputs,disinhibited-global_i, alpha=0.2, c='b', marker='.')
    ax3.set_ylabel('global feedback inhibition')
    ax3.set_xlabel('input R')
    ax3.set_ylim([-0.02,0.23])

    ax4.scatter(inputs,disinhibited, alpha=0.2, c='b', marker='.')
    ax4.set_ylabel('disinhibited')
    ax4.set_xlabel('input R')    

def show_MeanSD():
    f2, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 2.5))
    f2.subplots_adjust(hspace=0.5)

    ax1.plot(x, ngr_mean, c='b')
    ax1.fill_between(x, ngr_mean-ngr_std, ngr_mean+ngr_std, color='b',
                     alpha=0.2)
    ax1.set_xlabel('pattern shift')

    ax2.plot(x, nnfr_mean, c='r')
    ax2.fill_between(x, nnfr_mean-nnfr_std, nnfr_mean+nnfr_std, color='r',
                     alpha=0.2)
    ax2.set_xlabel('pattern shift')

    ax3.plot(x, ntr_mean, c='k')
    ax3.fill_between(x, ntr_mean-ntr_std, ntr_mean+ntr_std, color='k',
                     alpha=0.2)
    ax3.set_xlabel('pattern shift')


def show_pairs():
    f3, (ax1, ax2) = plt.subplots(1, 2, figsize=(6, 3))
    f3.subplots_adjust(hspace=0.5)
    ax1.plot(x, ngr_mean, c='b')
    ax1.fill_between(x, ngr_mean-ngr_std, ngr_mean+ngr_std, color='b',
                     alpha=0.2)
    ax1.plot(x, ntr_mean, c='k')
    ax1.fill_between(x, ntr_mean-ntr_std, ntr_mean+ntr_std, color='k',
                     alpha=0.2)

    ax2.plot(x, nnfr_mean, c='r')
    ax2.fill_between(x, nnfr_mean-nnfr_std, nnfr_mean+nnfr_std, color='r',
                     alpha=0.2)
    ax2.plot(x, ntr_mean, c='k')
    ax2.fill_between(x, ntr_mean-ntr_std, ntr_mean+ntr_std, color='k',
                     alpha=0.2)

show_scatter()