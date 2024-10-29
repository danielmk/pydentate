"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
from pydentate import net_basket_cell_ring, neuron_tools, oscillations_analysis
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


dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, 'output', 'figure_3_input_strength')
    
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
     # 'shuffling_coherence': []
     'area_over_line': [],
     'area_over_line_normalized': [],
     'input_strength': [],
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
    d['input_strength'].append(curr_data.root.parameters.input_current.read())
    
    """Calculate Powers"""
    freq = curr_data.root.analysis.fft.freq.read()
    amp = curr_data.root.analysis.fft.amp.read()
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
    
    # n_shuffles=10
    # shuffle = oscillations_analysis.shuffling_coherence(curr_data.root.times.read(), curr_data.root.parameters.dt.read(), curr_data.root.parameters.duration.read(), n_shuffles=n_shuffles)
    # d['shuffling_coherence'].append(np.abs(shuffle[0] - shuffle[1]).sum())

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

"""PLOTTING"""
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 22})

min_synapses_idx = df['input_strength'].argmin()
max_synapses_idx = df['input_strength'].argmax()

min_synapses_file = df['file'][min_synapses_idx]
max_synapses_file = df['file'][max_synapses_idx]
def exp_decay(t, tau, V):
    '''Exponential decay helper function'''
    return V * np.e ** (-t / tau)

"""Create the synaptic Kernel"""
decay_tau = 0.0018
t = np.arange(0,0.020, 0.0001)
kernel = exp_decay(t, decay_tau, 1)
kernel = kernel / kernel.sum()

curr_dt = 0.1
curr_binary = spike_to_x.units_times_to_binary(np.array(max_synapses_file.root.units.read(), dtype=int), max_synapses_file.root.times.read(),
                                               n_units=curr_data.root.parameters.n_pvbcs.read(),
                                               dt=curr_dt,
                                               total_time=curr_data.root.parameters.duration.read() * 1000 + 1)

curr_binary_convolved = np.array([np.convolve(x, kernel, 'same') for x in curr_binary])

synaptic_field = curr_binary_convolved.sum(axis=0)

t = np.arange(0, 2001.0, 0.1)

t_conv_full = np.arange(0, min_synapses_file.root.analysis.coherence.synaptic_field.read().shape[0] * 0.1, 0.1)

fig, ax = plt.subplots(3, 3)
ax[0, 0].scatter(min_synapses_file.root.times, min_synapses_file.root.units, marker='|', color=colors[3])
ax[0, 0].set_xlim((0, 200))
ax[0, 0].set_xlabel("Time (ms)")
ax[0, 0].set_ylabel("Neuron #")
ax[0, 0].set_title("Mean $f_{\mu}$:\n" + f"{min_synapses_file.root.analysis.frequencies.read().mean():.2f}")

ax[0, 1].scatter(max_synapses_file.root.times, max_synapses_file.root.units, marker='|', color=colors[4])
ax[0, 1].set_xlim((0, 200))
ax[0, 1].set_xlabel("Time (ms)")
ax[0, 1].set_ylabel("Neuron #")
ax[0, 1].set_title("Mean $f_{\mu}$:\n" + f"{max_synapses_file.root.analysis.frequencies.read().mean():.2f}")

ax[1, 0].plot(t_conv_full, min_synapses_file.root.analysis.coherence.synaptic_field.read(), color=colors[3])
ax[1, 0].set_xlim((0, 200))
#ax[1, 0].set_ylim((0, 2))
ax[1, 0].set_xlabel("Time (ms)")
ax[1, 0].set_ylabel("Synaptic Field")

ax[1, 1].plot(t_conv_full, max_synapses_file.root.analysis.coherence.synaptic_field.read(), color=colors[4])
ax[1, 1].set_xlim((0, 200))
#ax[1, 1].set_ylim((0, 2))
ax[1, 1].set_xlabel("Time (ms)")
ax[1, 1].set_ylabel("Synaptic Field")

ax[2, 0].plot(bin_size_ms, min_synapses_file.root.analysis.coherence.mean_k.read(), color=colors[3], marker='o')
ax[2, 1].plot(bin_size_ms, max_synapses_file.root.analysis.coherence.mean_k.read(), color=colors[4], marker='o')
ax[2, 0].set_xlabel("Bin $\\tau$ (ms)")
ax[2, 0].set_ylabel("Coherence Measure $k(\\tau)$")
ax[2, 1].set_xlabel("Bin $\\tau$ (ms)")
ax[2, 1].set_ylabel("Coherence Measure $k(\\tau)$")
dt = min_synapses_file.root.parameters.dt.read()
ax[2, 1].set_xlim((0, 8))
ax[2, 1].set_ylim((0, 1.3))

async_binary = spike_to_x.units_times_to_binary(
    np.array(min_synapses_file.root.units.read(), dtype=int),
    min_synapses_file.root.times.read(),
    n_units=min_synapses_file.root.parameters.n_pvbcs.read(),
    dt=dt,
    total_time=min_synapses_file.root.parameters.duration.read() * 1000 + 1)

sync_binary = spike_to_x.units_times_to_binary(
    np.array(max_synapses_file.root.units.read(), dtype=int),
    max_synapses_file.root.times.read(),
    n_units=max_synapses_file.root.parameters.n_pvbcs.read(),
    dt=dt,
    total_time=max_synapses_file.root.parameters.duration.read() * 1000 + 1)

async_peak_tau = (1 / df.iloc[min_synapses_idx]['mean_frequency']) * 1000  # in ms
sync_peak_tau = (1 / df.iloc[max_synapses_idx]['mean_frequency']) * 1000  # in ms

async_peak_coherence = np.nanmean(oscillations_analysis.pairwise_coherence(async_binary, bin_size=int(async_peak_tau / dt)))

sync_peak_coherence = np.nanmean(oscillations_analysis.pairwise_coherence(sync_binary, bin_size=int(sync_peak_tau / dt)))

ax[2, 0].plot([0, async_peak_tau], [0, async_peak_coherence], color='k', marker='o', linestyle='dashed')

ax[2, 1].plot([0, sync_peak_tau], [0, sync_peak_coherence], color='k', marker='o', linestyle='dashed')

sns.scatterplot(data=df, x='input_strength', y='mean_pearson_r', ax=ax[0, 2], color='k')

sns.scatterplot(data=df, x='input_strength', y='k_point_1_over_f', ax=ax[1, 2], color='k')

sns.scatterplot(data=df, x='input_strength', y='area_over_line', ax=ax[2, 2], color='k')


# linearity_n_points = 30
# oscillations_analysis.linearity_analysis(curr_binary, curr_dt, curr_duration, n_points=linearity_n_points)

fig, ax = plt.subplots(1, 4)
sns.scatterplot(data=df, x='input_strength', y='small_delta_coherence', ax=ax[0], color='k')
sns.scatterplot(data=df, x='input_strength', y='one_ms_delta_coherence', ax=ax[1], color='k')
sns.scatterplot(data=df, x='input_strength', y='k_point_1_over_f', ax=ax[2], color='k')
sns.scatterplot(data=df, x='input_strength', y='one_ms_delta_slope_coherence', ax=ax[3], color='k')

"""
plt.figure()
sns.boxplot(data=df, x='inhibition_on', y='k_point_1_over_f', hue='gap_junctions')

plt.figure()
sns.boxplot(data=df, x='inhibition_on', y='total_coherence', hue='gap_junctions')

plt.figure()
sns.boxplot(data=df, x='inhibition_on', y='k_1_over_f', hue='gap_junctions')
"""

"""
example_file = "pvring_seed_pp-weight_rec-weight_input-rate_n-pvbcs_n-input-syns_n-inputs_gap-resistance_gap-junctions_145_0.001_0.0076_2_120_100_192_600.0_True.h5"

example_path = os.path.join(results_dir, example_file)

file = tables.open_file(example_path)

times = file.root.times.read()
units = file.root.units.read()
duration = file.root.parameters.duration.read()
dt = file.root.parameters.dt.read()

n_timepoints = int((duration * 1000) / dt)

spike_counts = np.unique(times, return_counts=True)

spike_occurences_indices = np.array(spike_counts[0] / dt, dtype=int)

spike_counts_full = np.zeros(n_timepoints+1)

spike_counts_full[spike_occurences_indices] = spike_counts[1]

gaussian_sigma = 1 / dt
gaussian_width = gaussian_sigma * 10
gaussian_kernel = windows.gaussian(gaussian_width, gaussian_sigma)
gaussian_kernel = gaussian_kernel / gaussian_kernel.sum()

spike_counts_full_convolved = convolve(spike_counts_full, gaussian_kernel)

# hist, bins = np.histogram(output_spiketimes, bins=1000)

sampling_rate = 1 / (dt / 1000)
# sampling_rate = 5000

sp = fftshift(fft(spike_counts_full_convolved))
freq = fftshift(fftfreq(len(spike_counts_full_convolved), 1/sampling_rate))

amp = np.abs(sp)

fig, ax = plt.subplots(2, 1)
ax[1].plot(freq[len(freq)//2:], amp[len(amp)//2:])
ax[1].set_title('FFT')
ax[1].set_xlabel('Frequency (Hz)')
ax[1].set_ylabel('Amplitude')
ax[1].set_xlim(0,300)
ax[1].set_ylim(0, 10000)

fig.tight_layout()

ax[0].scatter(curr_times, curr_units, marker='|')
"""
