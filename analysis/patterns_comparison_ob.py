# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from pathlib2 import Path
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
basePath = Path(r"Z:\pyDentate\pyDentateData")
model_input = r"pattern_separation_data_local_input_revised"
seed = "seed10000"
scale = "scale1000"
filename = r"1_NDP_matrix.txt"


case = r"input_patterns"
file_to_open = basePath / model_input / seed / case / filename
data = np.loadtxt(str(file_to_open))
linear_data= np.array([])
for x in range(25):
    linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
inputs = linear_data

case = r"net_disinhibitedrev"
file_to_open = basePath / model_input / seed / scale /case / filename
try:
    data = np.loadtxt(str(file_to_open))
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
    data = np.loadtxt(str(file_to_open))
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
    data = np.loadtxt(str(file_to_open))
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
    data = np.loadtxt(str(file_to_open))
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
    data = np.loadtxt(str(file_to_open))
    linear_data= np.array([])
    for x in range(25):
        linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
    nonfacilitating = linear_data
except:
    nonfacilitating = tuned
    print("no nonfacilitating data")
    
case = r"net_reshuffledrev"
file_to_open = basePath / model_input / seed / scale /case / filename
try:
    data = np.loadtxt(str(file_to_open))
    linear_data= np.array([])
    for x in range(25):
        linear_data = np.concatenate((linear_data, data[x,x:data.shape[0]]))
    reshuffled = linear_data
except:
    reshuffled = tuned
    print("no reshuffled data")    

def show_scatter():
    f1, ((ax1, ax2, ax3, ax4, ax5), (ax6, ax7, ax8, ax9, ax10)) = plt.subplots(2, 5, figsize=(11, 5))
    f1.subplots_adjust(hspace=0.7, wspace=0.4)

    comparison1 = inputs - tuned
    ax1.scatter(inputs, comparison1, alpha=0.2, c='b', marker='.')
    ax1.set_title('full network', fontsize=10)
    ax1.set_ylabel('delta output R')
    ax1.set_xlabel('input corr.')
    ax1.set_ylim([-0.1, 0.5])
    significant, marks, means, bins = bindata(inputs, comparison1)
    ax1.scatter(significant, marks*0.45, c='k', marker='.')
    ax1.scatter(bins, means, c='k', marker='_')
    ax1.plot([0,1],[0,0],'k:')

    comparison2 = inputs - noFeedBack
    ax2.scatter(inputs, comparison2, alpha=0.2, c='b', marker='.')
    ax2.set_title('only FF', fontsize=10)
    ax2.set_xlabel('input R')
    ax2.set_ylim([-0.1, 0.5])
    significant, marks, means, bins = bindata(inputs, comparison2)
    ax2.scatter(significant, marks*0.45, c='k', marker='.')
    ax2.scatter(bins, means, c='k', marker='_')
    ax2.plot([0, 1], [0, 0], 'k:')

    comparison3 = inputs - disinhibited
    ax3.scatter(inputs, comparison3, alpha=0.2, c='b', marker='.')
    ax3.set_title('disinh.', fontsize=10)
    ax3.set_xlabel('input R')
    ax3.set_ylim([-0.1, 0.5])
    significant, marks, means, bins = bindata(inputs, comparison3)
    ax3.scatter(significant, marks*0.45, c='k', marker='.')
    ax3.scatter(bins, means, c='k', marker='_')
    ax3.plot([0, 1], [0, 0], 'k:')

    x_bar = [1, 2, 3]
    bar1 = comparison1[(inputs < 1) & (inputs > 0)]
    bar2 = comparison2[(inputs < 1) & (inputs > 0)]
    bar3 = comparison3[(inputs < 1) & (inputs > 0)]
    mean_diff = [np.mean(bar1), np.mean(bar2), np.mean(bar3)]
    SDs = [np.std(bar1), np.std(bar2), np.std(bar3)]
    SEMs = SDs/np.sqrt(325)
    ax4.bar(x_bar, mean_diff, yerr=SEMs, color='b', alpha=0.7)
    p1 = round(wilcoxon(bar1)[1]*3, 4)
    p2 = round(wilcoxon(bar2)[1]*3, 4)
    p3 = round(wilcoxon(bar3)[1]*3, 4)
    ax4.text(0.4, 0.2, "p=" + str(p1) + " p=" + str(p2) + " p=" + str(p3))

    comparison4 = tuned-global_i
    ax6.scatter(inputs, comparison4, alpha=0.2, c='b', marker='.')
    ax6.set_title('global feedback inh.', fontsize=10)
    ax6.set_xlabel('input R')
    ax6.set_ylim([-0.5, 0.25])
    significant, marks, means, bins = bindata(inputs, comparison4)
    ax6.scatter(significant, marks*0.2, c='k', marker='.')
    ax6.scatter(bins, means, c='k', marker='_')
    ax6.plot([0, 1], [0, 0], 'k:')

    comparison5 = tuned - nonfacilitating
    ax7.scatter(inputs, comparison5, alpha=0.2, c='b', marker='.')
    ax7.set_title('nonfaciliting', fontsize=10)
    ax7.set_xlabel('input R')
    ax7.set_ylim([-0.5, 0.25])
    significant, marks, means, bins = bindata(inputs, comparison5)
    ax7.scatter(significant, marks*0.2, c='k', marker='.')
    ax7.scatter(bins, means, c='k', marker='_')
    ax7.plot([0, 1], [0, 0], 'k:')

    comparison6 = tuned - reshuffled
    ax8.scatter(inputs, comparison6, alpha=0.2, c='b', marker='.')
    ax8.set_title('reshuffled', fontsize=10)
    ax8.set_xlabel('input R')
    ax8.set_ylim([-0.5, 0.25])
    significant, marks, means, bins = bindata(inputs, comparison6)
    ax8.scatter(significant, marks*0.2, c='k', marker='.')
    ax8.scatter(bins, means, c='k', marker='_')
    ax8.plot([0, 1], [0, 0], 'k:')

    x_bar = [1, 2, 3]
    bar1 = comparison4[(inputs < 1) & (inputs > 0)]
    bar2 = comparison5[(inputs < 1) & (inputs > 0)]
    bar3 = comparison6[(inputs < 1) & (inputs > 0)]
    mean_diff = [np.mean(bar1), np.mean(bar2), np.mean(bar3)]
    SDs = [np.std(bar1), np.std(bar2), np.std(bar3)]
    SEMs = SDs/np.sqrt(325)
    ax9.bar(x_bar, mean_diff, yerr=SEMs, color='b', alpha=0.7)
    p1 = round(wilcoxon(bar1)[1]*3, 4)
    p2 = round(wilcoxon(bar2)[1]*3, 4)
    p3 = round(wilcoxon(bar3)[1]*3, 4)
    ax9.text(0.4, 0.003, "p=" + str(p1) + " p=" + str(p2) + " p=" + str(p3))

    savename = 'comparison_' + model_input + seed + scale + '.pdf'
    f1.savefig(savename)

def bindata(inputs, linear_data):
    data = linear_data
    bins = np.linspace(0.05, 0.95, 10)
    digitized = np.digitize(inputs, bins)
    binned_data = [data[digitized == i] for i in range(len(bins))]
    p = [wilcoxon(binned_data[i][:])[1] for i in range(len(bins))]
    means = np.array([data[digitized == i].mean() for i in range(len(bins))])
    SDs = [data[digitized == i].std() for i in range(len(bins))]
    
    x = bins
    # plt.plot(x, means, color='b')
    # plt.fill_between(x, means - SDs, means + SDs, alpha=0.2, color='b')
    significant = x[np.array(p)<0.05/10]

    marks = np.ones(len(significant))
    return significant, marks, means, bins
    # plt.scatter(significant, marks, color='k', marker='.')

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