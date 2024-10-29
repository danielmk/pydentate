"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
from pydentate import net_basket_cell_ring, neuron_tools
from pydentate.inputs import inhom_poiss, inhomogeneous_poisson_process, sigmoid
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


dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, 'output', 'figure_5_slow_gamma_input_high_sync')
    
all_files = [os.path.join(data_dir, x) for x in os.listdir(data_dir)]

bin_size_ms = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

neuron_types = {'gc': 10000, 'bc': 120, 'hc': 120, 'mc': 300}

neuron_dict = {}

input_dict = {}

conv_dict = {}

def exp_decay(t, tau, V):
    '''Exponential decay helper function'''
    return V * np.e ** (-t / tau)

decay_tau = 0.0018
t = np.arange(0,0.020, 0.0001)
kernel = exp_decay(t, decay_tau, 1)
kernel = kernel / kernel.sum()

for cell in neuron_types.keys():
    d = {'seed': [],
         'rec_weight': [],
         'input_rate': [],
         'n_pvbcs': [],
         'n_input_syns': [],
         'n_inputs': [],
         'gap_resistance': [],
         'gap_junctions': [],
         'mean_frequency': [],
         'theta_power': [],  # 4-12 Hz
         'gamma_power': [],  # 20-100 Hz
         'low_gamma_power': [],  # 20-40
         'high_gamma_power': [], # 60-100
         'ripple_power': [],  # 120-250
         'total_coherence': [],
         'small_delta_coherence': [],
         'intermediate_delta_coherence': [],
         'file': [],
         'mean_rec_synapses': [],
         'k_point_1_over_f': [],
         'k_1_over_f': [],
         'one_ms_delta_coherence': [],
         'one_ms_delta_slope_coherence': [],
         'area_over_line': [],
         'area_over_line_normalized': [],
         'mean_k': [],
         'mean_pearson_r': []
         }

    for curr_file in all_files:
        curr_data = tables.open_file(curr_file, mode='r')
        fname = curr_file.split(os.path.sep)[-1]
        seed = curr_data.root.parameters.seed.read()
    
        d['seed'].append(curr_data.root.parameters.seed.read())
        d['rec_weight'].append(curr_data.root.parameters.rec_weight.read())
        d['input_rate'].append(curr_data.root.parameters.input_rate.read())
        d['n_pvbcs'].append(curr_data.root.parameters.n_pvbcs.read())
        d['n_input_syns'].append(curr_data.root.parameters.n_input_syns.read())
        d['n_inputs'].append(curr_data.root.parameters.n_inputs.read())
        d['gap_resistance'].append(curr_data.root.parameters.gap_resistance.read())
        d['gap_junctions'].append(curr_data.root.parameters.gap_junction.read())
        d['mean_frequency'].append(eval(f'curr_data.root.analysis.frequencies_{cell}.read().mean()'))
        d['mean_rec_synapses'].append(curr_data.root.parameters.n_rec_synapses.read().mean())
        d['mean_k'].append(eval(f'curr_data.root.analysis.coherence.mean_k_{cell}.read()'))
        
        """Calculate Powers"""
        freq = eval(f'curr_data.root.analysis.fft.freq_{cell}.read()')
        amp = eval(f'curr_data.root.analysis.fft.amp_{cell}.read()')
        amp = amp / amp.sum()
        theta_band = np.logical_and(12 > freq, freq > 4)
        gamma_band = np.logical_and(100 > freq, freq > 20)
        low_gamma_band = np.logical_and(40 > freq, freq > 20)
        high_gamma_band = np.logical_and(100 > freq, freq > 60)
        ripple_band = np.logical_and(250 > freq, freq > 120)
    
        d['theta_power'].append(amp[theta_band].mean())
        d['gamma_power'].append(amp[gamma_band].mean())
        d['low_gamma_power'].append(amp[low_gamma_band].mean())
        d['high_gamma_power'].append(amp[high_gamma_band].mean())
        d['ripple_power'].append(amp[ripple_band].mean())
        
        d['total_coherence'].append(np.nanmean(np.array(eval(f'curr_data.root.analysis.coherence.mean_k_{cell}.read()'))))
        d['small_delta_coherence'].append(eval(f'curr_data.root.analysis.coherence.mean_k_{cell}.read()[0]'))
        d['one_ms_delta_coherence'].append(eval(f'curr_data.root.analysis.coherence.mean_k_{cell}.read()[1]'))
        d['one_ms_delta_slope_coherence'].append(eval(f'curr_data.root.analysis.coherence.mean_k_{cell}.read()[2] - curr_data.root.analysis.coherence.mean_k_{cell}.read()[1]'))
        d['intermediate_delta_coherence'].append(eval(f'curr_data.root.analysis.coherence.mean_k_{cell}.read()[5]'))
        d['file'].append(curr_data)
    
        d['k_point_1_over_f'].append(float(eval(f'curr_data.root.analysis.coherence.k_point_1_over_f_{cell}.read()')))
        d['k_1_over_f'].append(float(eval(f'curr_data.root.analysis.coherence.k_1_over_f_{cell}.read()')))
    
        d['area_over_line'].append(float(eval(f'curr_data.root.analysis.coherence.area_over_line_{cell}.read()')))
        d['area_over_line_normalized'].append(float(eval(f'curr_data.root.analysis.coherence.area_over_line_normalized_{cell}.read()')))
        
        d['mean_pearson_r'].append(float(eval(F'curr_data.root.analysis.coherence.mean_pearson_r_{cell}.read()')))
        
        seed = curr_data.root.parameters.seed.read()
        
        np.random.seed(seed)

        duration = 2.0  # In s
        input_rate = 15
        n_inputs = 120
        peak_freq = 40
        dt = 0.1

        bins = np.arange(0, 2000.1, 0.1)

        temporal_patterns = np.array([inhomogeneous_poisson_process(0.0, duration, dt/1000, 30, peak_freq, refractory_period=0.001) for x in range(n_inputs)], dtype=object) * 1000
        
        temp_flat = np.hstack(temporal_patterns)
        
        temporal_binned, _ = np.histogram(temp_flat, bins)
        
        temporal_binned_convolved = np.convolve(temporal_binned, kernel, 'full')
        
        conv_dict[seed] = temporal_binned_convolved
        
        input_dict[seed] = temporal_patterns
        
        neuron_dict[f'{cell}'] = d.copy()
        
        """
        curr_data.create_array('/analysis/coherence', 'delta_t_point_1', obj=delta_t_point_1)
        curr_data.create_array('/analysis/coherence', 'k_point_1_over_f', obj=k_point_one_over_f.mean())
        
        delta_t_1 = (1 / frequencies.mean()) * 1000
        k_one_over_f = oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(delta_t_1 / curr_dt))
        
        curr_data.create_array('/analysis/coherence', 'delta_t_1', obj=delta_t_1)
        curr_data.create_array('/analysis/coherence', 'k_1_over_f', obj=k_one_over_f.mean())
        """


for cell in neuron_types.keys():
    df = pd.DataFrame(neuron_dict[cell])
    df['inhibition_on'] = df['rec_weight'] > 0
    
    new_category = []
    for i in df.index:
        if df.loc[i]['gap_junctions'] == True:
            string = 'on'
        else:
            string = 'off'
        if df.loc[i]['inhibition_on'] == True:
            string = string + '_on'
        else:
            string = string + '_off'
        new_category.append(string)
    new_category
    df['network_condition'] = pd.Series(new_category, dtype="category")
    
    """PLOTTING"""
    colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']
    plt.rcParams['svg.fonttype'] = 'none'
    font = {'family' : 'Arial',
            'size'   : 22}
    matplotlib.rc('font', **font)
    """
    fig, ax = plt.subplots(1, 4)
    sns.boxplot(data=df, x='network_condition', y='small_delta_coherence', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], ax=ax[0], palette=colors, dodge=False)
    sns.boxplot(data=df, x='network_condition', y='one_ms_delta_coherence', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=ax[1], palette=colors, dodge=False)
    sns.boxplot(data=df, x='network_condition', y='k_point_1_over_f', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=ax[2], palette=colors, dodge=False)
    sns.boxplot(data=df, x='network_condition', y='one_ms_delta_slope_coherence', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=ax[3], palette=colors, dodge=False)
    
    # gj, syn
    off_off = (df['gap_junctions'] == False) & (df['inhibition_on'] == False)
    off_on = (df['gap_junctions'] == False) & (df['inhibition_on'] == True)
    on_off = (df['gap_junctions'] == True) & (df['inhibition_on'] == False)
    on_on = (df['gap_junctions'] == True) & (df['inhibition_on'] == True)
    
    off_off_mean_k = np.array([x for x in df[off_off]['mean_k']])
    off_on_mean_k = np.array([x for x in df[off_on]['mean_k']])
    on_off_mean_k = np.array([x for x in df[on_off]['mean_k']])
    on_on_mean_k = np.array([x for x in df[on_on]['mean_k']])


    plt.figure()
    plt.errorbar(x=bin_size_ms, y=on_off_mean_k.mean(axis=0), yerr=off_off_mean_k.std(axis=0), marker='o', color=colors[0])
    plt.errorbar(x=bin_size_ms, y=on_off_mean_k.mean(axis=0), yerr=on_off_mean_k.std(axis=0), marker='o', color=colors[1])
    plt.errorbar(x=bin_size_ms, y=off_on_mean_k.mean(axis=0), yerr=off_on_mean_k.std(axis=0), marker='o', color=colors[2])
    plt.errorbar(x=bin_size_ms, y=on_on_mean_k.mean(axis=0), yerr=on_on_mean_k.std(axis=0), marker='o', color=colors[3])
    plt.legend(("off_off", "on_off", "off_on", "on_on"))
    plt.xlabel("Bin Tau (ms)")
    plt.ylabel("Coherence K(tau)")
    """
    
    off_off = (df['gap_junctions'] == False) & (df['inhibition_on'] == False)
    off_on = (df['gap_junctions'] == False) & (df['inhibition_on'] == True)
    on_off = (df['gap_junctions'] == True) & (df['inhibition_on'] == False)
    on_on = (df['gap_junctions'] == True) & (df['inhibition_on'] == True)
    
    off_off_mean_k = np.array([x for x in df[off_off]['mean_k']])
    off_on_mean_k = np.array([x for x in df[off_on]['mean_k']])
    on_off_mean_k = np.array([x for x in df[on_off]['mean_k']])
    on_on_mean_k = np.array([x for x in df[on_on]['mean_k']])

    # Plot Examples
    # gj, syn
    off_off_example = 2
    off_on_example = 0
    on_off_example = 3
    on_on_example = 1
    
    fig, ax = plt.subplots(3, 5)
    fig.suptitle(f"Celltype: {cell}")
    ax[0, 0].eventplot(input_dict[df.loc[off_off_example]['seed']], color='k')
    ax[0, 1].scatter(eval(f"df.loc[off_off_example]['file'].root.times.times_{cell}[df.loc[off_off_example]['file'].root.times.times_{cell}.read() < 2000]"), eval(f"df.loc[off_off_example]['file'].root.units.units_{cell}[df.loc[off_off_example]['file'].root.times.times_{cell}.read() < 2000]"), marker='|', color=colors[0])
    ax[0, 2].scatter(eval(f"df.loc[on_off_example]['file'].root.times.times_{cell}[df.loc[on_off_example]['file'].root.times.times_{cell}.read() < 2000]"), eval(f"df.loc[on_off_example]['file'].root.units.units_{cell}[df.loc[on_off_example]['file'].root.times.times_{cell}.read() < 2000]"), marker='|', color=colors[1])
    ax[0, 3].scatter(eval(f"df.loc[off_on_example]['file'].root.times.times_{cell}[df.loc[off_on_example]['file'].root.times.times_{cell}.read() < 2000]"), eval(f"df.loc[off_on_example]['file'].root.units.units_{cell}[df.loc[off_on_example]['file'].root.times.times_{cell}.read() < 2000]"), marker='|', color=colors[2])
    ax[0, 4].scatter(eval(f"df.loc[on_on_example]['file'].root.times.times_{cell}[df.loc[on_on_example]['file'].root.times.times_{cell}.read() < 2000]"), eval(f"df.loc[on_on_example]['file'].root.units.units_{cell}[df.loc[on_on_example]['file'].root.times.times_{cell}.read() < 2000]"), marker='|', color=colors[3])

    for a in ax[0, :]:
        a.set_xlabel("Time (ms)")

    ax[0,0].set_ylabel("Neuron #")

    for a in ax[0, :]:
        a.set_xlim((0, 2000))
    
    n_timepoints = df.loc[off_off_example]['file'].root.analysis.coherence.synaptic_field_gc.read().shape[0]
    t = np.arange(0, n_timepoints * 0.1, 0.1)
    t_conv = np.arange(0, conv_dict[df.loc[off_off_example]['seed']].shape[0] * 0.1, 0.1)
    
    ax[1, 0].plot(t_conv, conv_dict[df.loc[off_off_example]['seed']], color='k')
    ax[1, 1].plot(t, eval(f"df.loc[off_off_example]['file'].root.analysis.coherence.synaptic_field_{cell}.read()"), color=colors[0])
    ax[1, 2].plot(t, eval(f"df.loc[on_off_example]['file'].root.analysis.coherence.synaptic_field_{cell}.read()"), color=colors[1])
    ax[1, 3].plot(t, eval(f"df.loc[off_on_example]['file'].root.analysis.coherence.synaptic_field_{cell}.read()"), color=colors[2])
    ax[1, 4].plot(t, eval(f"df.loc[on_on_example]['file'].root.analysis.coherence.synaptic_field_{cell}.read()"), color=colors[3])
     
    for x in ax[1, :]:
        x.set_xlim((0, 2000))
        x.set_xlabel("Time (ms)")
        # x.set_ylim((0, 2.8))
        
    ax[1,0].set_ylabel("Synaptic Field")

    
    sns.boxplot(data=df, x='network_condition', y='mean_frequency', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], ax=ax[2, 0], palette=colors, dodge=False)
    
    """
    ax[1, 1].errorbar(x=bin_size_ms, y=off_off_mean_k.mean(axis=0), yerr=off_off_mean_k.std(axis=0), marker='o', color=colors[0])
    ax[1, 1].errorbar(x=bin_size_ms, y=on_off_mean_k.mean(axis=0), yerr=on_off_mean_k.std(axis=0), marker='o', color=colors[1])
    ax[1, 1].errorbar(x=bin_size_ms, y=off_on_mean_k.mean(axis=0), yerr=off_on_mean_k.std(axis=0), marker='o', color=colors[2])
    ax[1, 1].errorbar(x=bin_size_ms, y=on_on_mean_k.mean(axis=0), yerr=on_on_mean_k.std(axis=0), marker='o', color=colors[3])
    ax[1, 1].legend(("off_off", "on_off", "off_on", "on_on"))
    ax[1, 1].set_xlabel("Bin Tau (ms)")
    ax[1, 1].set_ylabel("Coherence K(tau)")
    """

    sns.boxplot(data=df, x='network_condition', y='k_point_1_over_f', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], ax=ax[2, 1], palette=colors, dodge=False)
    sns.boxplot(data=df, x='network_condition', y='area_over_line', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=ax[2, 2], palette=colors, dodge=False)
    sns.boxplot(data=df, x='network_condition', y='mean_pearson_r', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=ax[2, 3], palette=colors, dodge=False)

    plt.figure()
    plt.suptitle(f"{cell}")
    sns.boxplot(data=df, x='network_condition', y='low_gamma_power', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], palette=colors, dodge=False)

    """STATISTICS"""
    print(f"STATS FOR CELLTYPE: {cell}")
    formula = "mean_frequency ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df).fit()
    table = anova_lm(point_one_ols)
    print(f"Mean Frequency Stats: \n {table}")
    
    formula = "k_point_1_over_f ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df).fit()
    table = anova_lm(point_one_ols)
    print(f'k 0.1/f stats: \n {table}')
    
    formula = "area_over_line ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df).fit()
    table = anova_lm(point_one_ols)
    print(f'Area over line stats: \n {table}')
    
    formula = "mean_pearson_r ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df).fit()
    table = anova_lm(point_one_ols)
    print(f'Mean Pearson R: \n {table}')
    
    formula = "low_gamma_power ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df).fit()
    table = anova_lm(point_one_ols)
    print(f'Low Gamma Power: \n {table}')
    
    formula = "high_gamma_power ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df).fit()
    table = anova_lm(point_one_ols)
    print(f'High Gamma Power: \n {table}')
    
    formula = "gamma_power ~ C(gap_junctions) * C(inhibition_on)"
    point_one_ols = ols(formula, data=df).fit()
    table = anova_lm(point_one_ols)
    print(f'Overall Gamma Power: \n {table}')






