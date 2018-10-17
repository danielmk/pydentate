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
model_input = r"pattern_separation_data_local_30Hz_input"
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
                 r'net_reshuffledrev']:
 
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
    f1, ax_array = plt.subplots(6, 7, figsize=(11, 10))
    f1.subplots_adjust(hspace=0.5, wspace=0.4)
    
    results = {}
    summary_bar_mean = {}
    summary_bar_CoV = {}
    for i, seed in enumerate(seeds):
        inputs = dictionary[seed]['inputs']
        tuned = dictionary[seed][r'net_tunedrev']
        global_i = dictionary[seed][r'net_globalrev']
        noFB = dictionary[seed][r'net_nofeedbackrev']
        nonfac = dictionary[seed][r'net_nonfacilitatingrev']
        reshuff = dictionary[seed][r'net_reshuffledrev']
        results[seed] = {}
        
        ax1 = ax_array[0, i]
        comparison1 = noFB -  tuned
        ax1.set_title(seed, fontsize=10)
        if i == 0:
            ax1.set_ylabel('isolated FB')
        ax1.set_ylim([-0.1, 0.25])
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, comparison1)
        ax1.plot(bins, means, color='b')
        ax1.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        # ax1.scatter(significant, marks*0.2, c='k', marker='.')
        # ax6.plot(bins, means, c='k', linestyle='--')
        ax1.plot([0, 1], [0, 0], 'k:')
        results[seed]['bins'] = bins
        results[seed]['FB'] = means
    
        ax2 = ax_array[1, i]
        norm = 1 # np.mean(comparison1[(inputs < 1) & (inputs > 0)])
        comparison2 = (global_i - tuned)/norm
        # ax7.scatter(inputs, comparison5, alpha=0.2, c='b', marker='.')
        if i == 0:
            ax2.set_ylabel('tuning effect', fontsize=10)
        # ax2.set_xlabel('input R')
        ax2.set_ylim([-0.1, 0.1])
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, comparison2)
        ax2.plot(bins, means, color='b')
        ax2.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        # ax2.scatter(significant, marks*0.1, c='k', marker='.')
        # ax7.plot(bins, means, c='k', linestyle='--')
        ax2.plot([0, 1], [0, 0], 'k:')
        results[seed]['tuning'] = means
    
        ax3 = ax_array[2, i]
        comparison3 = (nonfac - tuned)/norm
        #ax8.scatter(inputs, comparison6, alpha=0.2, c='b', marker='.')
        if i == 0:
            ax3.set_ylabel('facilitation effect', fontsize=10)
        # ax3.set_xlabel('input R')
        ax3.set_ylim([-0.1, 0.1])
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, comparison3)
        ax3.plot(bins, means, color='b')
        ax3.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        # ax3.scatter(significant, marks*0.1, c='k', marker='.')
        #ax8.plot(bins, means, c='k', linestyle='--')
        ax3.plot([0, 1], [0, 0], 'k:')
        results[seed]['facilitation'] = means
        
        ax4 = ax_array[3, i]
        comparison4 = (reshuff - tuned)/norm
        if i == 0:
            ax4.set_ylabel('shuffling effect', fontsize=10)
        ax4.set_xlabel('input R')
        ax4.set_ylim([-0.1, 0.1])
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, comparison4)
        ax4.plot(bins, means, color='b')
        ax4.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax4.plot([0, 1], [0, 0], 'k:')
        results[seed]['shuffling'] = means

        ax5 = ax_array[4, i]
        x_bar = [1, 2, 3]
        # bar1 = comparison2[(inputs < 1) & (inputs > 0)]
        # bar2 = comparison3[(inputs < 1) & (inputs > 0)]
        # bar3 = comparison4[(inputs < 1) & (inputs > 0)]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, comparison2)
        bar1 = means[1:-1]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, comparison3)
        bar2 = means[1:-1]
        significant, marks, means, SDs, SEMs, bins = bindata(inputs, comparison4)
        bar3 = means[1:-1]
        mean_diff = [np.mean(bar1), np.mean(bar2), np.mean(bar3)]
        SDs = [np.std(bar1), np.std(bar2), np.std(bar3)]
        SEMs = SDs/np.sqrt(len(comparison2))
        ax5.bar(x_bar, mean_diff, yerr=SEMs, color='b', alpha=0.7)
        ax5.set_xticklabels(['','tuning','facilitation','reshuffling'], rotation= 30, fontsize=8)
        summary_bar_mean[seed] = mean_diff
        
        ax6 = ax_array[5, i]
        significant, marks, means, tunedSD, SEMs, bins = bindata(inputs, tuned)
        tunedCoV = tunedSD/(bins-means)
        ind = [0,1,2,9,10,11]   # these values are removed from CoV because they are highly unrliable due to small mean
        tunedCoV[ind] = np.nan
        ax6.plot(bins, tunedCoV, c='b')
        results[seed]['tunedCoV'] = tunedCoV 
        ax6.set_ylim([0, 0.5])
        significant, marks, means, globalSD, SEMs, bins = bindata(inputs, global_i)
        globalCoV = globalSD/(bins-means)
        globalCoV[ind] = np.nan
        ax6.plot(bins, globalCoV, c='k')
        results[seed]['globalCoV'] = globalCoV

        significant, marks, means, nonfacSD, SEMs, bins = bindata(inputs, nonfac)
        nonfacCoV = nonfacSD/(bins-means)
        nonfacCoV[ind] = np.nan
        results[seed]['nonfacCoV'] = nonfacCoV 
        significant, marks, means, reshuffSD, SEMs, bins = bindata(inputs, reshuff)
        reshuffCoV = reshuffSD/(bins-means)
        reshuffCoV[ind] = np.nan
        results[seed]['reshuffCoV'] = reshuffCoV
        CoV = [1 - np.mean(globalCoV[3:-3])/np.mean(tunedCoV[3:-3]),
               1 - np.mean(nonfacCoV[3:-3])/np.mean(tunedCoV[3:-3]),
               1 - np.mean(reshuffCoV[3:-3])/np.mean(tunedCoV[3:-3])]
        summary_bar_CoV[seed] = CoV
        
    savename = 'robustness_' + model_input + scale + '.pdf'
    f1.savefig(savename)
    
    f2, ((ax1, ax2, ax3, ax4, ax5), 
         (ax6, ax7, ax8, ax9, ax10)) = plt.subplots(2, 5, figsize=(11, 4))
    f2.subplots_adjust(hspace=0.5, wspace=0.4)
    p = pd.DataFrame.from_dict(results, orient='index')
    
    ax1.set_title('isolated FB', fontsize=10)
    ax1.set_ylim([-0.01, 0.15])
    # ax1.set_xlabel('input R')
    bins = p.bins.mean()
    means = p.FB.mean()
    SEMs = p.FB.apply(pd.Series).std()/np.sqrt(7)
    ax1.plot(bins, means, color='b')
    ax1.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
    ax1.plot([0, 1], [0, 0], 'k:')

    ax2.set_title('tuning effect', fontsize=10)
    # ax2.set_xlabel('input R')
    ax2.set_ylim([-0.05, 0.05])
    means = p.tuning.mean()
    SEMs = p.tuning.apply(pd.Series).std()/np.sqrt(7)
    ax2.plot(bins, means, color='b')
    ax2.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
    # ax2.scatter(significant, marks*0.1, c='k', marker='.')
    ax2.plot([0, 1], [0, 0], 'k:')

    ax3.set_title('facilitation effect', fontsize=10)
    # ax3.set_xlabel('input R')
    ax3.set_ylim([-0.05, 0.05])
    means = p.facilitation.mean()
    SEMs = p.facilitation.apply(pd.Series).std()/np.sqrt(7)
    ax3.plot(bins, means, color='b')
    ax3.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
    # ax3.scatter(significant, marks*0.1, c='k', marker='.')
    ax3.plot([0, 1], [0, 0], 'k:')

    ax4.set_title('shuffling effect', fontsize=10)
    # ax4.set_xlabel('input R')
    ax4.set_ylim([-0.05, 0.05])
    means = p.shuffling.mean()
    SEMs = p.shuffling.apply(pd.Series).std()/np.sqrt(7)
    ax4.plot(bins, means, color='b')
    ax4.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
    # ax4.scatter(significant, marks*0.1, c='k', marker='.')
    ax4.plot([0, 1], [0, 0], 'k:')

    x_bar = [1, 2, 3]
    means = pd.DataFrame.from_dict(summary_bar_mean, orient='index')
    SEMs = pd.DataFrame.from_dict(summary_bar_mean, orient='index').std()/np.sqrt(7)
    mean_diff = means.mean()
    ax5.bar(x_bar, mean_diff, yerr=SEMs, color='b', alpha=0.7)
    ax5.set_xticklabels(['','tuning','facilitation','reshuffling'], rotation= 30)
    ax5.set_ylim([-0.05, 0.05])

    savetxtname = 'grandmean_' + model_input + scale + '.txt'
    np.savetxt(savetxtname, means)
    
    # CoVs
    ax6.set_title('tuned vs global CoV', fontsize=10)
    #ax6.set_xlabel('input R')
    ax6.set_ylabel('CoV')
    gmeans = p.globalCoV.mean()
    tmeans = p.tunedCoV.mean()
    ax6.plot(bins, tmeans, color='b')
    ax6.plot(bins, gmeans,  color='k')
    gSEMs = p.globalCoV.apply(pd.Series).std()/np.sqrt(7)
    tSEMs = p.tunedCoV.apply(pd.Series).std()/np.sqrt(7)
    ax6.fill_between(bins, tmeans - tSEMs, tmeans + tSEMs, alpha=0.2, color='b')
    ax6.fill_between(bins, gmeans - gSEMs, gmeans + gSEMs, alpha=0.2, color='k')
    ax6.plot([0, 1], [0, 0], 'k:')
    # ax6.set_ylim([-0.01, 0.1])
    
    ax7.set_title('tuning CoV effect', fontsize=10)
    # ax7.set_xlabel('input R')
    # ax7.set_ylabel('CoV')
    diff_CoV = p.tunedCoV - p.globalCoV
    diffCoVmeans = diff_CoV.mean()
    ax7.plot(bins, diffCoVmeans, color='b')
    diffCoV_SEM = diff_CoV.apply(pd.Series).std()/np.sqrt(7)
    ax7.fill_between(bins, diffCoVmeans - diffCoV_SEM,
                     diffCoVmeans + diffCoV_SEM, alpha=0.2, color='b')
    ax7.set_ylim([-0.1, 0.1])    
    ax7.plot([0, 1], [0, 0], 'k:')

    ax8.set_title('facilitation CoV effect', fontsize=10)
    #ax8.set_xlabel('input R')
    # ax8.set_ylabel('CoV')
    diff_CoV = p.tunedCoV - p.nonfacCoV
    diffCoVmeans = diff_CoV.mean()
    ax8.plot(bins, diffCoVmeans, color='b')
    diffCoV_SEM = diff_CoV.apply(pd.Series).std()/np.sqrt(7)
    ax8.fill_between(bins, diffCoVmeans - diffCoV_SEM,
                     diffCoVmeans + diffCoV_SEM, alpha=0.2, color='b')
    ax8.set_ylim([-0.1, 0.1])    
    ax8.plot([0, 1], [0, 0], 'k:')
    
    ax9.set_title('shuffling CoV effect', fontsize=10)
    #ax6.set_xlabel('input R')
    # ax9.set_ylabel('CoV')
    # ax6.set_ylim([0, 0.15])
    diff_CoV = p.tunedCoV - p.reshuffCoV
    diffCoVmeans = diff_CoV.mean()
    ax9.plot(bins, diffCoVmeans, color='b')
    diffCoV_SEM = diff_CoV.apply(pd.Series).std()/np.sqrt(7)
    ax9.fill_between(bins, diffCoVmeans - diffCoV_SEM,
                     diffCoVmeans + diffCoV_SEM, alpha=0.2, color='b')
    ax9.set_ylim([-0.1, 0.1])    
    ax9.plot([0, 1], [0, 0], 'k:')

    x_bar = [1, 2, 3]
    means = pd.DataFrame.from_dict(summary_bar_CoV, orient='index')
    SEMs = pd.DataFrame.from_dict(summary_bar_CoV, orient='index').std()/np.sqrt(7)
    CoV_diff = means.mean()
    ax10.bar(x_bar, CoV_diff, yerr=SEMs, color='b', alpha=0.7)
    ax10.set_xticklabels(['','tuning','facilitation','reshuffling'], rotation= 30)
    ax10.set_ylim([-0.2, 0.2])
    
    savename = 'gand_mean_' + model_input + scale + '.pdf'
    savetxtname = 'CoV_' + model_input + scale + '.txt'
    np.savetxt(savetxtname, means)
    f2.savefig(savename)


def bindata(inputs, linear_data):
    data = linear_data
    bins = np.linspace(0, 1, 11)
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

show_scatter()
