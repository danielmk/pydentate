# -*- coding: utf-8 -*-
"""
Spyder Editor

time resolved analysis for seeds 10000 to 10007.
data points are the mean values within each bin of each seed
(10 bins excluding 0 and 1 for input R; 8 seeds)
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
# from scipy.stats import binned_statistic
from scipy.stats import wilcoxon
import matplotlib.colors as colors
# import pandas as pd
# from scipy.optimize import curve_fit

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# what to analyse
basePath = Path(r"R:\pyDentate\pyDentateData")
model_input = r"pattern_separation_data_local_input_revised"
scale = "scale1000"
filename = r"1_leutgeb-measure-tresolved_len-bin_1000.txt"

savename = 'comparison_thetaResolved' + model_input + scale + '.pdf'

dictionary = {}
seeds = ["seed10000", "seed10001", "seed10002", "seed10003", "seed10004",
         "seed10005", "seed10006"]

timebins = 5

for seed in seeds:
    dictionary[seed] = {}
    for case in [r"input_patterns", r"net_disinhibitedrev",
                 r"net_nofeedbackrev", r"net_tunedrev", r"net_globalrev",
                 r"net_nonfacilitatingrev", r"net_reshuffledrev"]:
        if case == r"input_patterns":
            file_to_open = basePath / model_input / seed / case / filename
        else:
            file_to_open = basePath / model_input / seed / scale /case / filename
        try:
            data = np.loadtxt(file_to_open)
            mask = np.array(range(timebins, np.shape(data)[1]))
            data = np.delete(data, mask, axis=1)
            dictionary[seed][case] = data
        except:
            dictionary[seed][case] = np.ones((325, timebins))
            print("no " + seed + " " + case)
    
def show():
    f1, ax_array = plt.subplots(5, timebins + 1, figsize=(timebins, 5))
    f1.subplots_adjust(hspace=0.5, wspace=0.4)
    
    full_effect = np.zeros((len(seeds), timebins, 2))
    FB_effect = np.zeros((len(seeds), timebins, 2))
    FF_effect = np.zeros((len(seeds), timebins, 2))
    tuning_effect = np.zeros((len(seeds), timebins, 2))
    facilitation_effect = np.zeros((len(seeds), timebins, 2))
    
    full_effect_im = np.zeros((11, timebins))
    FB_effect_im = np.zeros((11, timebins))
    FF_effect_im = np.zeros((11, timebins))
    tuning_effect_im = np.zeros((11, timebins))
    facilitation_effect_im = np.zeros((11, timebins))
    
    full_effect_per_cycle = np.zeros((11, timebins, len(seeds)))
    FB_effect_per_cycle = np.zeros((11, timebins, len(seeds)))
    
    inputs = np.zeros([len(seeds), len(data), timebins])
    tuned = np.zeros([len(seeds), len(data), timebins])
    noFB = np.zeros([len(seeds), len(data), timebins])
    disinh = np.zeros([len(seeds), len(data), timebins])
    glob = np.zeros([len(seeds), len(data), timebins])
    nonfac = np.zeros([len(seeds), len(data), timebins])
    
    for i, seed in enumerate(seeds):
        inputs[i,:,:] = (dictionary[seed][r"input_patterns"])
        if np.isnan(inputs).any():
            print(seed)
            print(np.count_nonzero(np.isnan(inputs)))
        np.nan_to_num(inputs)
        tuned[i,:,:] = (dictionary[seed][r"net_tunedrev"])
        noFB[i,:,:] = (dictionary[seed][r"net_nofeedbackrev"])
        disinh[i,:,:] = (dictionary[seed][r"net_disinhibitedrev"])
        glob[i,:,:] = (dictionary[seed][r"net_globalrev"])
        nonfac[i,:,:] = (dictionary[seed][r"net_nonfacilitatingrev"])
        
    
    for t in range(timebins):
        # full PS effect
        ax1 = ax_array[0, t]
        c1 = np.zeros([len(seeds), 12])
        x = np.zeros([len(seeds), 12])
        for s in range(len(seeds)):
            inp = inputs[s, :, t]
            comparison1 = inp - tuned[s, :, t]
            significant, marks, means, SDs, SEMs, bins = bindata(inp, comparison1)
            c1[s] = means
            x[s] = bins
            # ax1.plot(bins, means, color='b', alpha=0.2)
            full_effect[s, t, 0] = np.nanmean(comparison1[(inp > 0) & (inp < 1)])
            # full_effect[s, t, 1] = np.nanstd(comparison1[(inp > 0) & (inp < 1)])
        # ax1.scatter(x, c1, alpha=0.2, c='b', marker='.')
        significant, marks, means, SDs, SEMs, bins = bindata(x, c1)
        full_effect_im[:, t] = means[1:]
        ax1.plot(bins, means, color='b')
        ax1.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax1.set_title('time bin ' + str(t), fontsize=10)
        if t == 0:
            ax1.set_ylabel('full PS effect')
        else:
            ax1.tick_params(labelleft='off')  
        ax1.set_ylim([-0.1, 0.7])
        # ax1.plot(x, c1, c='k', linestyle='--')
        ax1.plot([0, 1], [0, 0], 'k:')
    
        # FB effect
        ax2 = ax_array[1, t]
        c2 = np.zeros([len(seeds), 12])
        for s in range(len(seeds)):
            inp = inputs[s, :, t]
            comparison2 = noFB[s, :, t] - tuned[s, :, t]
            significant, marks, means, SDs, SEMs, bins = bindata(inp, comparison2)
            # ax2.plot(bins, means, alpha=0.2, c='b')
            c2[s] = means
            FB_effect_per_cycle[:, t, s] = means[1:]
            FB_effect[s, t, 0] = np.nanmean(comparison2[(inp > 0) & (inp < 1)])
            FB_effect[s, t, 1] = np.nanstd(comparison2[(inp > 0) & (inp < 1)])     
        if t == 0:
            ax2.set_ylabel('FB effect', fontsize=10)
        else:
            ax2.tick_params(labelleft='off')
        significant, marks, means, SDs, SEMs, bins = bindata(x, c2)
        FB_effect_im[:, t] = means[1:]
        ax2.plot(bins, means, color='b')
        ax2.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax2.set_ylim([-0.1, 0.7])
        # ax2.plot(bins, means, c='k', linestyle='--')
        ax2.plot([0, 1], [0, 0], 'k:')
        
    
        # FF effect
        ax3 = ax_array[2, t]
        c3 = np.zeros([len(seeds), 12])
        for s in range(len(seeds)):
            inp = inputs[s, :, t]
            comparison3 = disinh[s, :, t] - noFB[s, :, t]
            significant, marks, means, SDs, SEMs, bins = bindata(inp, comparison3)
            c3[s] = means
            FF_effect[s, t, 0] = np.nanmean(comparison3[(inp > 0) & (inp < 1)])
            FF_effect[s, t, 1] = np.nanstd(comparison3[(inp > 0) & (inp < 1)])
            # ax3.plot(bins, means, alpha=0.2, c='b')
        if t == 0:
            ax3.set_ylabel('FF effect', fontsize=10)
        else:
            ax3.tick_params(labelleft='off')  
        significant, marks, means, SDs, SEMs, bins = bindata(x, c3)
        FF_effect_im[:, t] = means[1:]
        ax3.plot(bins, means, color='b')
        ax3.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax3.set_ylim([-0.1, 0.2])
        # ax3.plot(bins, means, c='k', linestyle='--')
        ax3.plot([0, 1], [0, 0], 'k:')
       
        ax4 = ax_array[3, t]
        c4 = np.zeros([len(seeds), 12])
        for s in range(len(seeds)):
            inp = inputs[s, :, t]
            norm = FB_effect[s, t, 0]
            comparison4 = glob[s, :, t] - tuned[s, :, t]
            significant, marks, means, SDs, SEMs, bins = bindata(inp, comparison4)
            c4[s] = means
            tuning_effect[s, t, 0] = np.nanmean(comparison4[(inp > 0) & (inp < 1)])
            tuning_effect[s, t, 1] = np.nanstd(comparison4[(inp > 0) & (inp < 1)])
            # ax4.plot(bins, means, alpha=0.2, c='b')
        if t == 0:
            ax4.set_ylabel('tuning effect', fontsize=10)
        else:
            ax4.tick_params(labelleft='off')  
        significant, marks, means, SDs, SEMs, bins = bindata(x, c4)
        tuning_effect_im[:, t] = means[1:]
        ax4.plot(bins, means, color='b')
        ax4.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax4.set_ylim([-0.1, 0.2])
        # ax4.plot(bins, means, color='b')
        # ax4.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax4.plot([0, 1], [0, 0], 'k:')
    
        ax5 = ax_array[4, t]
        c5 = np.zeros([len(seeds), 12])
        for s in range(len(seeds)):
            inp = inputs[s, :, t]
            norm = FB_effect[s, t, 0]
            comparison5 = nonfac[s, :, t] - tuned[s, :, t]
            significant, marks, means, SDs, SEMs, bins = bindata(inp, comparison5)
            c5[s] = means
            facilitation_effect[s, t, 0] = np.nanmean(comparison5[(inp > 0) & (inp < 1)])
            facilitation_effect[s, t, 1] = np.nanstd(comparison5[(inp > 0) & (inp < 1)])
            # ax5.plot(bins, means, alpha=0.2, c='b')
        if t == 0:
            ax5.set_ylabel('facilitation effect', fontsize=10)
        else:
            ax5.tick_params(labelleft='off')  
        significant, marks, means, SDs, SEMs, bins = bindata(x, c5)
        facilitation_effect_im[:, t] = means[1:]
        ax5.plot(bins, means, color='b')
        ax5.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax5.set_ylim([-0.1, 0.2])
        # ax5.set_ylim([-0.2, 0.2])
        # ax5.plot(bins, means, color='b')
        # ax5.fill_between(bins, means - SEMs, means + SEMs, alpha=0.2, color='b')
        ax5.plot([0, 1], [0, 0], 'k:')
        
    
    x_bar = np.array(range(timebins))
    ax = ax_array[0, timebins]
    means = full_effect[:, :, 0].mean(axis=0)
    SEMs = full_effect[:,:,0].std(axis=0)/np.sqrt(7)
    ax.plot(x_bar, means, color='b')
    ax.fill_between(x_bar, means - SEMs, means + SEMs, alpha=0.2, color='b')
    ax.plot([0, timebins], [0, 0], 'k:')
    ax.set_ylim([-0.1, 0.7])

    ax = ax_array[1, timebins]
    means = FB_effect[:, :, 0].mean(axis=0)
    SEMs = FB_effect[:,:,0].std(axis=0)/np.sqrt(7)
    ax.plot(x_bar, means, color='b')
    ax.fill_between(x_bar, means - SEMs, means + SEMs, alpha=0.2, color='b')
    ax.plot([0, timebins], [0, 0], 'k:')
    ax.set_ylim([-0.1, 0.7])

    ax = ax_array[2, timebins]
    means = FF_effect[:, :, 0].mean(axis=0)
    SEMs = FF_effect[:,:,0].std(axis=0)/np.sqrt(7)
    ax.plot(x_bar, means, color='b')
    ax.fill_between(x_bar, means - SEMs, means + SEMs, alpha=0.2, color='b')
    ax.plot([0, timebins], [0, 0], 'k:')
    ax.set_ylim([-0.1, 0.2])

    ax = ax_array[3, timebins]
    means = tuning_effect[:, :, 0].mean(axis=0)
    SEMs = tuning_effect[:,:,0].std(axis=0)/np.sqrt(7)
    ax.plot(x_bar, means, color='b')
    ax.fill_between(x_bar, means - SEMs, means + SEMs, alpha=0.2, color='b')
    ax.plot([0, timebins], [0, 0], 'k:')
    ax.set_ylim([-0.1, 0.2])

    ax = ax_array[4, timebins]
    means = facilitation_effect[:, :, 0].mean(axis=0)
    SEMs = facilitation_effect[:,:,0].std(axis=0)/np.sqrt(7)
    ax.plot(x_bar, means, color='b')
    ax.fill_between(x_bar, means - SEMs, means + SEMs, alpha=0.2, color='b')
    ax.plot([0, timebins], [0, 0], 'k:')
    ax.set_ylim([-0.1, 0.2])

    f1.savefig(savename)
    
    f2, ax_array = plt.subplots(1, 5, figsize=(timebins, 2))
    f2.subplots_adjust(hspace=0.5, wspace=0.4)
    
    # cNorm = colors.Normalize(vmin=np.min(full_effect_im), vmax=np.max(full_effect_im))
    cNorm = colors.Normalize(vmin=0, vmax=0.2)
    cmap = 'jet'
    
    ax = ax_array[0]
    ax.imshow(full_effect_im, cmap=cmap, norm=cNorm)
    ax.set_title('full effect')
    ax.set_ylabel('input R', fontsize=10)
    ax.set_xlabel('cycle')
    
    ax = ax_array[1]
    ax.imshow(FB_effect_im, cmap=cmap, norm=cNorm)
    ax.set_title('FB effect')
    ax.set_xlabel('cycle')
    
    ax = ax_array[2]
    ax.imshow(FF_effect_im, cmap=cmap, norm=cNorm)
    ax.set_title('FF effect')
    ax.set_xlabel('cycle')
    
    ax = ax_array[3]
    ax.imshow(tuning_effect_im, cmap=cmap, norm=cNorm)
    ax.set_title('tuning effect')
    ax.set_xlabel('cycle')
    
    ax = ax_array[4]
    ax.imshow(facilitation_effect_im, cmap=cmap, norm=cNorm)
    ax.set_title('facilitation effect')
    ax.set_xlabel('cycle')
    
    cbar_ax = f2.add_axes([0.92, 0.15, 0.01, 0.7])
    mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=cNorm,
                              orientation='vertical', label='PS effect')
    name = 'heatmap_thetaresolved' + model_input + '.pdf'
    f2.savefig(name)
    name = 'thetaresolved_fulleffect_' + model_input + '.csv'
    np.savetxt(name, full_effect_im)
    name = 'thetaresolved_FBeffect_' + model_input + '.csv'
    np.savetxt(name, FB_effect_im)

def bindata(inputs, linear_data):
    data = linear_data
    bins = np.linspace(0, 1, 11)
    digitized = np.digitize(inputs, bins)
    binned_data = [data[digitized == i] for i in range(len(bins))]
    # p = [wilcoxon(binned_data[i][:])[1] for i in range(len(bins))]
    means = np.zeros(len(bins)+1)
    SDs = np.zeros(len(bins)+1)
    Ns = np.zeros(len(bins)+1)
    for i in range(len(bins)+1):
        means[i] = np.nanmean(data[digitized == i])
        SDs[i] = np.nanstd(data[digitized == i])
        Ns[i] = sum(~np.isnan(data[digitized == i]))
    # means = np.array([data[digitized == i].mean() for i in range(len(bins))])
    # SDs = [data[digitized == i].std() for i in range(len(bins))]
    # n = [len(binned_data[i][:]) for i in range(len(bins))]
    SEMs = SDs/np.sqrt(Ns)
    # plt.plot(x, means, color='b')
    # plt.fill_between(x, means - SDs, means + SDs, alpha=0.2, color='b')
        # plt.scatter(significant, marks, color='k', marker='.')
    # significant = bins[np.array(p) < 0.05/10]
    # marks = np.ones(len(significant))
    significant = []
    marks = []
    # means = np.append(0, means)
    # SDs = np.append(0, SDs)
    # SEMs = np.append(0, SEMs)
    bins = bins + 0.05  # bin middles for plotting
    bins[-1] = 1
    bins = np.append(0, bins)

    return significant, marks, means, SDs, SEMs, bins

show()