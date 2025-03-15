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
import matplotlib as mpl



dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, 'output', 'figure_2_on_off')
    
all_files = [os.path.join(data_dir, x) for x in os.listdir(data_dir)]

bin_size_ms = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

d = {'seed': [],
     'pp_weight': [],
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
     'mean_pearson_r': []
     }

for curr_file in all_files:
    curr_data = tables.open_file(curr_file, mode='r')
    fname = curr_file.split(os.path.sep)[-1]

    d['seed'].append(curr_data.root.parameters.seed.read())
    d['pp_weight'].append(curr_data.root.parameters.pp_bc_weight.read())
    d['rec_weight'].append(curr_data.root.parameters.rec_weight.read())
    d['input_rate'].append(curr_data.root.parameters.input_rate.read())
    d['n_pvbcs'].append(curr_data.root.parameters.n_pvbcs.read())
    d['n_input_syns'].append(curr_data.root.parameters.n_input_syns.read())
    d['n_inputs'].append(curr_data.root.parameters.n_inputs.read())
    d['gap_resistance'].append(curr_data.root.parameters.gap_resistance.read())
    d['gap_junctions'].append(curr_data.root.parameters.gap_junction.read())
    d['mean_frequency'].append(curr_data.root.analysis.frequencies.read().mean())
    d['mean_rec_synapses'].append(curr_data.root.parameters.n_rec_synapses.read().mean())
    
    """Calculate Powers"""
    freq = curr_data.root.analysis.fft.freq.read()
    amp = curr_data.root.analysis.fft.amp.read()
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
    d['total_coherence'].append(np.nanmean(np.array(curr_data.root.analysis.coherence.mean_k.read())))
    d['small_delta_coherence'].append(curr_data.root.analysis.coherence.mean_k.read()[0])
    d['one_ms_delta_coherence'].append(curr_data.root.analysis.coherence.mean_k.read()[1])
    d['one_ms_delta_slope_coherence'].append(curr_data.root.analysis.coherence.mean_k.read()[2] - curr_data.root.analysis.coherence.mean_k.read()[1])
    d['intermediate_delta_coherence'].append(curr_data.root.analysis.coherence.mean_k.read()[5])
    d['file'].append(curr_data)

    d['k_point_1_over_f'].append(float(curr_data.root.analysis.coherence.k_point_1_over_f.read()))
    d['k_1_over_f'].append(float(curr_data.root.analysis.coherence.k_1_over_f.read()))

    d['area_over_line'].append(float(curr_data.root.analysis.coherence.area_over_line.read()))
    d['area_over_line_normalized'].append(float(curr_data.root.analysis.coherence.area_over_line_normalized.read()))
    d['mean_pearson_r'].append(float(curr_data.root.analysis.coherence.mean_pearson_r.read()))
    
    
    """
    curr_data.create_array('/analysis/coherence', 'delta_t_point_1', obj=delta_t_point_1)
    curr_data.create_array('/analysis/coherence', 'k_point_1_over_f', obj=k_point_one_over_f.mean())
    
    delta_t_1 = (1 / frequencies.mean()) * 1000
    k_one_over_f = oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(delta_t_1 / curr_dt))
    
    curr_data.create_array('/analysis/coherence', 'delta_t_1', obj=delta_t_1)
    curr_data.create_array('/analysis/coherence', 'k_1_over_f', obj=k_one_over_f.mean())
    """


df = pd.DataFrame(d)
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
mpl.rc('font', **font)

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

off_off_mean_k = np.array([x.root.analysis.coherence.mean_k.read() for x in df[off_off]['file']])
off_on_mean_k = np.array([x.root.analysis.coherence.mean_k.read() for x in df[off_on]['file']])
on_off_mean_k = np.array([x.root.analysis.coherence.mean_k.read() for x in df[on_off]['file']])
on_on_mean_k = np.array([x.root.analysis.coherence.mean_k.read() for x in df[on_on]['file']])

plt.figure()
plt.errorbar(x=bin_size_ms, y=off_off_mean_k.mean(axis=0), yerr=off_off_mean_k.std(axis=0), marker='o', color=colors[0])
plt.errorbar(x=bin_size_ms, y=on_off_mean_k.mean(axis=0), yerr=on_off_mean_k.std(axis=0), marker='o', color=colors[1])
plt.errorbar(x=bin_size_ms, y=off_on_mean_k.mean(axis=0), yerr=off_on_mean_k.std(axis=0), marker='o', color=colors[2])
plt.errorbar(x=bin_size_ms, y=on_on_mean_k.mean(axis=0), yerr=on_on_mean_k.std(axis=0), marker='o', color=colors[3])
plt.legend(("off_off", "on_off", "off_on", "on_on"))
plt.xlabel("Bin Tau (ms)")
plt.ylabel("Coherence K(tau)")

# Plot Examples
# gj, syn
off_off_example = 2
off_on_example = 0
on_off_example = 3
on_on_example = 1


fig, ax = plt.subplots(3, 4)
ax[0, 0].scatter(df.loc[off_off_example]['file'].root.times[df.loc[off_off_example]['file'].root.times.read() < 200], df.loc[off_off_example]['file'].root.units[df.loc[off_off_example]['file'].root.times.read() < 200], marker='|', color=colors[0])
ax[0, 1].scatter(df.loc[on_off_example]['file'].root.times[df.loc[on_off_example]['file'].root.times.read() < 200], df.loc[on_off_example]['file'].root.units[df.loc[on_off_example]['file'].root.times.read() < 200], marker='|', color=colors[1])
ax[0, 2].scatter(df.loc[off_on_example]['file'].root.times[df.loc[off_on_example]['file'].root.times.read() < 200], df.loc[off_on_example]['file'].root.units[df.loc[off_on_example]['file'].root.times.read() < 200], marker='|', color=colors[2])
ax[0, 3].scatter(df.loc[on_on_example]['file'].root.times[df.loc[on_on_example]['file'].root.times.read() < 200], df.loc[on_on_example]['file'].root.units[df.loc[on_on_example]['file'].root.times.read() < 200], marker='|', color=colors[3])

for a in ax[0, :]:
    a.set_xlim((0, 200))

ax[0,0].set_ylabel("Neuron #")


n_timepoints = df.loc[off_off_example]['file'].root.analysis.coherence.synaptic_field.read().shape[0]
t = np.arange(0, n_timepoints * 0.1, 0.1)

ax[1, 0].plot(t, df.loc[off_off_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[0])
ax[1, 1].plot(t, df.loc[on_off_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[1])
ax[1, 2].plot(t, df.loc[off_on_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[2])
ax[1, 3].plot(t, df.loc[on_on_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[3])

for x in ax[1, :]:
    x.set_xlim((0, 200))
    x.set_xlabel("Time (ms)")
    x.set_ylim((0, 2.8))
    
ax[1,0].set_ylabel("Synaptic Field")

# fig, ax = plt.subplots(1, 5)

sns.boxplot(data=df, x='network_condition', y='mean_frequency', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], ax=ax[2, 0], palette=colors, dodge=False)
"""
ax[1].errorbar(x=bin_size_ms, y=off_off_mean_k.mean(axis=0), yerr=off_off_mean_k.std(axis=0), marker='o', color=colors[0])
ax[1].errorbar(x=bin_size_ms, y=on_off_mean_k.mean(axis=0), yerr=on_off_mean_k.std(axis=0), marker='o', color=colors[1])
ax[1].errorbar(x=bin_size_ms, y=off_on_mean_k.mean(axis=0), yerr=off_on_mean_k.std(axis=0), marker='o', color=colors[2])
ax[1].errorbar(x=bin_size_ms, y=on_on_mean_k.mean(axis=0), yerr=on_on_mean_k.std(axis=0), marker='o', color=colors[3])
ax[1].legend(("off_off", "on_off", "off_on", "on_on"))
ax[1].set_xlabel("Bin Tau (ms)")
ax[1].set_ylabel("Coherence K(tau)")
"""

sns.boxplot(data=df, x='network_condition', y='k_point_1_over_f', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], ax=ax[2, 1], palette=colors, dodge=False)
sns.boxplot(data=df, x='network_condition', y='area_over_line', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=ax[2, 2], palette=colors, dodge=False)
sns.boxplot(data=df, x='network_condition', y='mean_pearson_r', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=ax[2, 3], palette=colors, dodge=False)

"""STATISTICS"""
formula = "mean_frequency ~ C(gap_junctions) * C(inhibition_on)"
point_one_ols = ols(formula, data=df).fit()
table = anova_lm(point_one_ols)
print(f"Mean Frequency Stats: \n {table.to_string()}")
# table.to_csv(r'C:\Users\Daniel\Dropbox\031_DG_Inhibition_Resonant_Properties\figures\Figure2\mean_frequency_stats.txt')

formula = "k_point_1_over_f ~ C(gap_junctions) * C(inhibition_on)"
point_one_ols = ols(formula, data=df).fit()
table = anova_lm(point_one_ols)
print(f'k 0.1/f stats: \n {table.to_string()}')
# table.to_csv(r'C:\Users\Daniel\Dropbox\031_DG_Inhibition_Resonant_Properties\figures\Figure2\k_point_1_over_f_stats.txt')

formula = "area_over_line ~ C(gap_junctions) * C(inhibition_on)"
point_one_ols = ols(formula, data=df).fit()
table = anova_lm(point_one_ols)
print(f'Area over line stats: \n {table.to_string()}')
# table.to_csv(r'C:\Users\Daniel\Dropbox\031_DG_Inhibition_Resonant_Properties\figures\Figure2\area_over_line_stats.txt')

formula = "mean_pearson_r ~ C(gap_junctions) * C(inhibition_on)"
point_one_ols = ols(formula, data=df).fit()
table = anova_lm(point_one_ols)
print(f'Mean Pearson R: \n {table.to_string()}')
# table.to_csv(r'C:\Users\Daniel\Dropbox\031_DG_Inhibition_Resonant_Properties\figures\Figure2\area_over_line_stats.txt')


# fig, ax = plt.subplots(3, 4)
fig = plt.figure()
spec = mpl.gridspec.GridSpec(4,4, figure=fig)
off_off_spikes = fig.add_subplot(spec[0, 0])
on_off_spikes = fig.add_subplot(spec[0, 1])
off_on_spikes = fig.add_subplot(spec[0, 2])
on_on_spikes = fig.add_subplot(spec[0, 3])

off_off_spikes.scatter(df.loc[off_off_example]['file'].root.times[df.loc[off_off_example]['file'].root.times.read() < 200], df.loc[off_off_example]['file'].root.units[df.loc[off_off_example]['file'].root.times.read() < 200], marker='|', color=colors[0])
on_off_spikes.scatter(df.loc[on_off_example]['file'].root.times[df.loc[on_off_example]['file'].root.times.read() < 200], df.loc[on_off_example]['file'].root.units[df.loc[on_off_example]['file'].root.times.read() < 200], marker='|', color=colors[1])
off_on_spikes.scatter(df.loc[off_on_example]['file'].root.times[df.loc[off_on_example]['file'].root.times.read() < 200], df.loc[off_on_example]['file'].root.units[df.loc[off_on_example]['file'].root.times.read() < 200], marker='|', color=colors[2])
on_on_spikes.scatter(df.loc[on_on_example]['file'].root.times[df.loc[on_on_example]['file'].root.times.read() < 200], df.loc[on_on_example]['file'].root.units[df.loc[on_on_example]['file'].root.times.read() < 200], marker='|', color=colors[3])

for a in [off_off_spikes, on_off_spikes, off_on_spikes, on_on_spikes]:
    a.set_xlim((0, 200))

off_off_spikes.set_ylabel("Neuron #")


n_timepoints = df.loc[off_off_example]['file'].root.analysis.coherence.synaptic_field.read().shape[0]
t = np.arange(0, n_timepoints * 0.1, 0.1)

off_off_field = fig.add_subplot(spec[1, 0])
on_off_field = fig.add_subplot(spec[1, 1])
off_on_field = fig.add_subplot(spec[1, 2])
on_on_field = fig.add_subplot(spec[1, 3])

off_off_field.plot(t, df.loc[off_off_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[0])
on_off_field.plot(t, df.loc[on_off_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[1])
off_on_field.plot(t, df.loc[off_on_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[2])
on_on_field.plot(t, df.loc[on_on_example]['file'].root.analysis.coherence.synaptic_field.read(), color=colors[3])

for x in [off_off_field, on_off_field, off_on_field, on_on_field]:
    x.set_xlim((0, 200))
    x.set_xlabel("Time (ms)")
    x.set_ylim((0, 2.8))
    
off_off_field.set_ylabel("Synaptic Field")


freq_plot = fig.add_subplot(spec[2,:2])
k_plot = fig.add_subplot(spec[2, 2:])
area_plot = fig.add_subplot(spec[3, :2])
pearson_plot = fig.add_subplot(spec[3, 2:])

sns.boxplot(data=df, x='network_condition', y='mean_frequency', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], ax=freq_plot, palette=colors, dodge=False)
sns.boxplot(data=df, x='network_condition', y='k_point_1_over_f', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'], ax=k_plot, palette=colors, dodge=False)
sns.boxplot(data=df, x='network_condition', y='area_over_line', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=area_plot, palette=colors, dodge=False)
sns.boxplot(data=df, x='network_condition', y='mean_pearson_r', hue='network_condition', order=['off_off', 'on_off', 'off_on', 'on_on'], hue_order=['off_off', 'on_off', 'off_on', 'on_on'],ax=pearson_plot, palette=colors, dodge=False)

"""CALCULATE VARIANCE AND STD"""
syn_stats = []

for f in d['file']:
    n_syn = f.root.parameters.chem_connection_matrix.read()
    mean = n_syn.sum(axis=0).mean()
    std = n_syn.sum(axis=0).std()
    var = n_syn.sum(axis=0).var()
    
    syn_stats.append((mean, std, var))
    
syn_stats = np.array(syn_stats)
idc_sorted = np.argsort(syn_stats[:, 0])

syn_stats = syn_stats[idc_sorted,:]
    
sd_color = 'r'
fig, ax1 = plt.subplots()
ax1.scatter(syn_stats[:,0], syn_stats[:,1], color='r')
ax1.set_xlabel("Mean N Syn")
ax1.set_ylabel("SD", color=sd_color)
ax1.tick_params(axis='y', labelcolor=sd_color)

color = 'tab:blue'
ax2 = ax1.twinx()
ax2.scatter(syn_stats[:,0], syn_stats[:,2], color=color)
ax2.set_ylabel("Variance", color=color)
ax2.tick_params(axis='y', labelcolor=color)