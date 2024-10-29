"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
from pydentate import net_basket_cell_ring, neuron_tools
from pydentate.inputs import inhom_poiss, homogeneous_poisson_process, sigmoid
import os
import argparse
import scipy.stats as stats
import platform
import matplotlib.pyplot as plt
from scipy.signal import windows, convolve
from scipy.ndimage import convolve
from scipy.spatial import distance_matrix
from numpy.fft import fft, fftshift, fftfreq
import sys
import tables
import networkx
from pydentate import spike_to_x
from pydentate import oscillations_analysis
import pdb
import pandas as pd
import seaborn as sns
from distutils.util import strtobool
from statsmodels.formula.api import ols
from statsmodels.graphics.api import interaction_plot, abline_plot
from statsmodels.stats.anova import anova_lm
import matplotlib
import copy


dirname = os.path.dirname(__file__)
homogeneous_dir = os.path.join(dirname, 'output', 'figure_6_full_dg_gc_membrane')
homogeneous_files = [os.path.join(homogeneous_dir, x) for x in os.listdir(homogeneous_dir)]

slow_gamma_dir = os.path.join(dirname, 'output', 'figure_6_slow_gamma_input_high_sync_gc_membrane_with_dendrite')
slow_gamma_files = [os.path.join(slow_gamma_dir, x) for x in os.listdir(slow_gamma_dir)]

fast_gamma_dir = os.path.join(dirname, 'output', 'figure_6_fast_gamma_input_high_sync_gc_membrane_with_dendrite')
fast_gamma_files = [os.path.join(fast_gamma_dir, x) for x in os.listdir(fast_gamma_dir)]

on_off_dict = {'off/off': [],
             'on/off': [],
             'off/on': [],
             'on/on': []
             }

freq_dict = {'homogeneous': copy.deepcopy(on_off_dict),
             'slow_gamma': copy.deepcopy(on_off_dict),
             'fast_gamma': copy.deepcopy(on_off_dict)}

sxx_dict = {'homogeneous': copy.deepcopy(on_off_dict),
             'slow_gamma': copy.deepcopy(on_off_dict),
             'fast_gamma': copy.deepcopy(on_off_dict)}

membrane_dict = {'homogeneous': copy.deepcopy(on_off_dict),
             'slow_gamma': copy.deepcopy(on_off_dict),
             'fast_gamma': copy.deepcopy(on_off_dict)}

theta_dict = {'homogeneous': copy.deepcopy(on_off_dict),
             'slow_gamma': copy.deepcopy(on_off_dict),
             'fast_gamma': copy.deepcopy(on_off_dict)}

gamma_dict = {'homogeneous': copy.deepcopy(on_off_dict),
             'slow_gamma': copy.deepcopy(on_off_dict),
             'fast_gamma': copy.deepcopy(on_off_dict)}

for input_type in ['homogeneous', 'slow_gamma', 'fast_gamma']:
    if input_type == 'homogeneous':
        all_files = homogeneous_files
    elif input_type == 'slow_gamma':
        all_files = slow_gamma_files
    elif input_type == 'fast_gamma':
        all_files = fast_gamma_files


    for curr_file in all_files:
        
        curr_data = tables.open_file(curr_file, mode='r')
        
        path_split = curr_file.split(os.path.sep)
        file_split = path_split[-1].split('_')
        
        if file_split[-2] == 'True':
            gj_on = True
        else:
            gj_on = False
        # gj_on = bool(file_split[-2])
        syn_on = bool(bool(float(file_split[-5])))
        
        print(f'gj: {file_split[-2]} {gj_on} syn: {file_split[-5]} {syn_on}')
        
        condition = (gj_on, syn_on)
        
        if condition == (False, False):
            freq = curr_data.root.frequency.read()

            freq_dict[input_type]['off/off'].append(curr_data.root.frequency.read())
            sxx_dict[input_type]['off/off'].append(curr_data.root.Pxx_den.read())
            if len(curr_data.root.gc_membrane.read().shape) > 1:
                membrane_dict[input_type]['off/off'].append(curr_data.root.gc_membrane.read().mean(axis=0))
            else:
                membrane_dict[input_type]['off/off'].append(curr_data.root.gc_membrane.read())
            
            theta_band = (freq > 4) & (freq < 12)
            
            gamma_band = (freq > 30) & (freq < 100)
            
            theta_dict[input_type]['off/off'].append(curr_data.root.Pxx_den.read()[theta_band].sum())
            
            gamma_dict[input_type]['off/off'].append(curr_data.root.Pxx_den.read()[gamma_band].sum())

        elif condition == (True, False):
            freq = curr_data.root.frequency.read()
            
            freq_dict[input_type]['on/off'].append(curr_data.root.frequency.read())
            sxx_dict[input_type]['on/off'].append(curr_data.root.Pxx_den.read())
            if len(curr_data.root.gc_membrane.read().shape) > 1:
                membrane_dict[input_type]['on/off'].append(curr_data.root.gc_membrane.read().mean(axis=0))
            else:
                membrane_dict[input_type]['on/off'].append(curr_data.root.gc_membrane.read())
            
            theta_band = (freq > 4) & (freq < 12)
            
            gamma_band = (freq > 30) & (freq < 100)
            
            theta_dict[input_type]['on/off'].append(curr_data.root.Pxx_den.read()[theta_band].sum())
            
            gamma_dict[input_type]['on/off'].append(curr_data.root.Pxx_den.read()[gamma_band].sum())

        elif condition == (False, True):
            freq = curr_data.root.frequency.read()
            
            freq_dict[input_type]['off/on'].append(curr_data.root.frequency.read())
            sxx_dict[input_type]['off/on'].append(curr_data.root.Pxx_den.read())
            if len(curr_data.root.gc_membrane.read().shape) > 1:
                membrane_dict[input_type]['off/on'].append(curr_data.root.gc_membrane.read().mean(axis=0))
            else:
                membrane_dict[input_type]['off/on'].append(curr_data.root.gc_membrane.read())
            
            theta_band = (freq > 4) & (freq < 12)
            
            gamma_band = (freq > 30) & (freq < 100)
            
            theta_dict[input_type]['off/on'].append(curr_data.root.Pxx_den.read()[theta_band].sum())
            
            gamma_dict[input_type]['off/on'].append(curr_data.root.Pxx_den.read()[gamma_band].sum())

        elif condition == (True, True):
            freq = curr_data.root.frequency.read()
            
            freq_dict[input_type]['on/on'].append(curr_data.root.frequency.read())
            sxx_dict[input_type]['on/on'].append(curr_data.root.Pxx_den.read())
            if len(curr_data.root.gc_membrane.read().shape) > 1:
                membrane_dict[input_type]['on/on'].append(curr_data.root.gc_membrane.read().mean(axis=0))
            else:
                membrane_dict[input_type]['on/on'].append(curr_data.root.gc_membrane.read())
            
            theta_band = (freq > 4) & (freq < 12)
            
            gamma_band = (freq > 30) & (freq < 100)
            
            theta_dict[input_type]['on/on'].append(curr_data.root.Pxx_den.read()[theta_band].sum())
            
            gamma_dict[input_type]['on/on'].append(curr_data.root.Pxx_den.read()[gamma_band].sum())


"""PLOTTING"""
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']
plt.rcParams['svg.fonttype'] = 'none'
font = {'family' : 'Arial',
        'size'   : 22}
matplotlib.rc('font', **font)

fig, ax = plt.subplots(4,3)

t = np.arange(0, 2000.2, 0.1)

for cond_idx, cond in enumerate(sxx_dict.keys()):
    for on_off_idx, on_off in enumerate(on_off_dict.keys()):
        ax[0, cond_idx].plot(t, membrane_dict[cond][on_off][0], color=colors[on_off_idx])
        ax[1, cond_idx].plot(freq_dict[cond][on_off][0], np.array(sxx_dict[cond][on_off]).mean(axis=0), color=colors[on_off_idx])
    ax[0, cond_idx].set_title(cond)
    ax[0, cond_idx].set_ylabel("Average GC Voltage (mV)")
    ax[0, cond_idx].set_xlabel("Time (ms)")
    # ax[0, cond_idx].legend(on_off_dict.keys())
    ax[1, cond_idx].set_xlim((0,100))
    # ax[1, cond_idx].legend(on_off_dict.keys())
    ax[1, cond_idx].set_xlabel("Frequency (Hz)")
    ax[1, cond_idx].set_ylabel("PSD (V$^2$/Hz)")
    
    df_theta = pd.DataFrame(theta_dict[cond])
    
    df_theta = df_theta.melt(var_name='condition')
    
    sns.boxplot(data=df_theta, x='condition', y='value', ax=ax[2, cond_idx], order=on_off_dict.keys(), hue_order=on_off_dict.keys(), palette=colors[:4])
    
    df_gamma = pd.DataFrame(gamma_dict[cond])
    
    df_gamma = df_gamma.melt(var_name='condition')
    
    sns.boxplot(data=df_gamma, x='condition', y='value', ax=ax[3, cond_idx], order=on_off_dict.keys(), hue_order=on_off_dict.keys(), palette=colors[:4])

    gj_on_theta = [x[0] for x in df_theta['condition'].astype(str).str.split('/')]
    
    inhibition_on_theta = [x[1] for x in df_theta['condition'].astype(str).str.split('/')]
    
    df_theta['gap_junctions'] = gj_on_theta
    
    df_theta['inhibition_on'] = inhibition_on_theta
    
    gj_on_gamma = [x[0] for x in df_gamma['condition'].astype(str).str.split('/')]
    
    inhibition_on_gamma = [x[1] for x in df_gamma['condition'].astype(str).str.split('/')]
    
    df_gamma['gap_junctions'] = gj_on_gamma
    
    df_gamma['inhibition_on'] = inhibition_on_gamma
    
    df_theta = df_theta.rename(columns={'value': 'theta_power'})
    
    df_gamma = df_gamma.rename(columns={'value': 'gamma_power'})
    print(f"Stats for {cond}\n")
    formula = "theta_power ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df_theta).fit()
    table = anova_lm(point_one_ols)
    print(f"Theta Stats: \n {table}")

    formula = "gamma_power ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df_gamma).fit()
    table = anova_lm(point_one_ols)
    print(f"Gamma Stats: \n {table}")
ax[0, 0].legend(on_off_dict.keys())



