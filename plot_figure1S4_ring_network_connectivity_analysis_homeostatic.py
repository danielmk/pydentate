"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import seaborn as sns
import tables
from scipy.signal import windows, convolve
from scipy.spatial import distance_matrix
from numpy.fft import fft, fftshift, fftfreq
from distutils.util import strtobool
from neuron import h, gui  # gui necessary for some parameters to h namespace
from pydentate import net_basket_cell_ring, neuron_tools, oscillations_analysis, spike_to_x
from pydentate.inputs import inhom_poiss, homogeneous_poisson_process, sigmoid


# Directories
dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, 'output', 'figure_1_S4_homeostatic_synapse_number')
all_files_homeostatic = [os.path.join(data_dir, x) for x in os.listdir(data_dir)]
data_dir = os.path.join(dirname, 'output', 'figure_1_connectivity')
all_files_normal = [os.path.join(data_dir, x) for x in os.listdir(data_dir)]

# Constants
BIN_SIZE_MS = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

# Initialize data dictionary
data_dict = {'normal':{
    'seed': [],
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
    'low_gamma_power': [],  # 20-40 Hz
    'high_gamma_power': [],  # 60-100 Hz
    'ripple_power': [],  # 120-250 Hz
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
    },
    'homeostatic':{
    'seed': [],
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
    'low_gamma_power': [],  # 20-40 Hz
    'high_gamma_power': [],  # 60-100 Hz
    'ripple_power': [],  # 120-250 Hz
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
}


# Process files
data_sources = {'normal': all_files_normal, 'homeostatic': all_files_homeostatic}
for data_source in data_sources:
    for curr_file in data_sources[data_source]:
        curr_data = tables.open_file(curr_file, mode='r')

        # Read parameters and analysis results
        params = curr_data.root.parameters
        analysis = curr_data.root.analysis
        coherence = analysis.coherence

        data_dict[data_source]['seed'].append(params.seed.read())
        data_dict[data_source]['pp_weight'].append(params.pp_bc_weight.read())
        data_dict[data_source]['rec_weight'].append(params.rec_weight.read())
        data_dict[data_source]['input_rate'].append(params.input_rate.read())
        data_dict[data_source]['n_pvbcs'].append(params.n_pvbcs.read())
        data_dict[data_source]['n_input_syns'].append(params.n_input_syns.read())
        data_dict[data_source]['n_inputs'].append(params.n_inputs.read())
        data_dict[data_source]['gap_resistance'].append(params.gap_resistance.read())
        data_dict[data_source]['gap_junctions'].append(params.gap_junction.read())
        data_dict[data_source]['mean_frequency'].append(analysis.frequencies.read().mean())
        data_dict[data_source]['mean_rec_synapses'].append(params.n_rec_synapses.read().mean())

        # Calculate Powers
        freq = analysis.fft.freq.read()
        amp = analysis.fft.amp.read() / analysis.fft.amp.read().sum()
        theta_band = np.logical_and(12 > freq, freq > 4)
        gamma_band = np.logical_and(100 > freq, freq > 20)
        low_gamma_band = np.logical_and(40 > freq, freq > 20)
        high_gamma_band = np.logical_and(100 > freq, freq > 60)
        ripple_band = np.logical_and(250 > freq, freq > 120)
    
        data_dict[data_source]['theta_power'].append(amp[theta_band].mean())
        data_dict[data_source]['gamma_power'].append(amp[gamma_band].mean())
        data_dict[data_source]['low_gamma_power'].append(amp[low_gamma_band].mean())
        data_dict[data_source]['high_gamma_power'].append(amp[high_gamma_band].mean())
        data_dict[data_source]['ripple_power'].append(amp[ripple_band].mean())
        data_dict[data_source]['total_coherence'].append(np.nanmean(np.array(coherence.mean_k.read())))
        data_dict[data_source]['small_delta_coherence'].append(coherence.mean_k.read()[0])
        data_dict[data_source]['one_ms_delta_coherence'].append(coherence.mean_k.read()[1])
        data_dict[data_source]['one_ms_delta_slope_coherence'].append(coherence.mean_k.read()[2] - coherence.mean_k.read()[1])
        data_dict[data_source]['intermediate_delta_coherence'].append(coherence.mean_k.read()[5])
        data_dict[data_source]['file'].append(curr_data)
        data_dict[data_source]['k_point_1_over_f'].append(float(coherence.k_point_1_over_f.read()))
        data_dict[data_source]['k_1_over_f'].append(float(coherence.k_1_over_f.read()))
        data_dict[data_source]['area_over_line'].append(float(coherence.area_over_line.read()))
        data_dict[data_source]['area_over_line_normalized'].append(float(coherence.area_over_line_normalized.read()))
        data_dict[data_source]['mean_pearson_r'].append(float(coherence.mean_pearson_r.read()))

# Create DataFrame
df = pd.DataFrame(data_dict['homeostatic'])
df['inhibition_on'] = df['rec_weight'] > 0

df_normal = pd.DataFrame(data_dict['normal'])

# Plotting
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 22})

# Get files with min and max synapses
min_synapses_idx = df['mean_rec_synapses'].argmin()
max_synapses_idx = df['mean_rec_synapses'].argmax()
min_synapses_file = df['file'][min_synapses_idx]
max_synapses_file = df['file'][max_synapses_idx]

def exp_decay(t, tau, V):
    """Exponential decay helper function"""
    return V * np.e ** (-t / tau)

# Create the synaptic kernel
decay_tau = 0.0018
t = np.arange(0, 0.020, 0.0001)
kernel = exp_decay(t, decay_tau, 1)
kernel = kernel / kernel.sum()

curr_dt = 0.1
curr_binary = spike_to_x.units_times_to_binary(
    np.array(max_synapses_file.root.units.read(), dtype=int),
    max_synapses_file.root.times.read(),
    n_units=max_synapses_file.root.parameters.n_pvbcs.read(),
    dt=curr_dt,
    total_time=max_synapses_file.root.parameters.duration.read() * 1000 + 1
)

curr_binary_convolved = np.array([np.convolve(x, kernel, 'same') for x in curr_binary])
synaptic_field = curr_binary_convolved.sum(axis=0)

t = np.arange(0, 2001.0, 0.1)
t_conv_full = np.arange(0, min_synapses_file.root.analysis.coherence.synaptic_field.read().shape[0] * 0.1, 0.1)

fig, ax = plt.subplots(3, 3)
ax[0, 0].scatter(min_synapses_file.root.times, min_synapses_file.root.units, marker='|', color=colors[3])
ax[0, 0].set_xlim((0, 200))
ax[0, 0].set_xlabel("Time (ms)")
ax[0, 0].set_ylabel("Neuron #")
ax[0, 0].set_title("Mean $N_{syn}$:\n" + f"{min_synapses_file.root.parameters.n_rec_synapses.read().mean():.2f}")

ax[0, 1].scatter(max_synapses_file.root.times, max_synapses_file.root.units, marker='|', color=colors[4])
ax[0, 1].set_xlim((0, 200))
ax[0, 1].set_xlabel("Time (ms)")
ax[0, 1].set_ylabel("Neuron #")
ax[0, 1].set_title("Mean $N_{syn}$:\n" + f"{max_synapses_file.root.parameters.n_rec_synapses.read().mean():.2f}")

ax[1, 0].plot(t_conv_full, min_synapses_file.root.analysis.coherence.synaptic_field.read(), color=colors[3])
ax[1, 0].set_xlim((0, 200))
ax[1, 0].set_ylim((0, 2))
ax[1, 0].set_xlabel("Time (ms)")
ax[1, 0].set_ylabel("Synaptic Field")

ax[1, 1].plot(t_conv_full, max_synapses_file.root.analysis.coherence.synaptic_field.read(), color=colors[4])
ax[1, 1].set_xlim((0, 200))
ax[1, 1].set_ylim((0, 2))
ax[1, 1].set_xlabel("Time (ms)")
ax[1, 1].set_ylabel("Synaptic Field")

ax[2, 0].plot(BIN_SIZE_MS, min_synapses_file.root.analysis.coherence.mean_k.read(), color=colors[3], marker='o')
ax[2, 1].plot(BIN_SIZE_MS, max_synapses_file.root.analysis.coherence.mean_k.read(), color=colors[4], marker='o')
ax[2, 0].set_xlabel("Bin $\\tau$ (ms)")
ax[2, 0].set_ylabel("Coherence Measure $k(\\tau)$")
ax[2, 1].set_xlabel("Bin $\\tau$ (ms)")
ax[2, 1].set_ylabel("Coherence Measure $k(\\tau)$")

dt = min_synapses_file.root.parameters.dt.read()

async_binary = spike_to_x.units_times_to_binary(
    np.array(min_synapses_file.root.units.read(), dtype=int),
    min_synapses_file.root.times.read(),
    n_units=min_synapses_file.root.parameters.n_pvbcs.read(),
    dt=dt,
    total_time=min_synapses_file.root.parameters.duration.read() * 1000 + 1
)

sync_binary = spike_to_x.units_times_to_binary(
    np.array(max_synapses_file.root.units.read(), dtype=int),
    max_synapses_file.root.times.read(),
    n_units=max_synapses_file.root.parameters.n_pvbcs.read(),
    dt=dt,
    total_time=max_synapses_file.root.parameters.duration.read() * 1000 + 1
)

async_peak_tau = (1 / df.iloc[min_synapses_idx]['mean_frequency']) * 1000  # in ms
sync_peak_tau = (1 / df.iloc[max_synapses_idx]['mean_frequency']) * 1000  # in ms

async_peak_coherence = np.nanmean(oscillations_analysis.pairwise_coherence(async_binary, bin_size=int(async_peak_tau / dt)))
sync_peak_coherence = np.nanmean(oscillations_analysis.pairwise_coherence(sync_binary, bin_size=int(sync_peak_tau / dt)))

ax[2, 0].plot([0, async_peak_tau], [0, async_peak_coherence], color='k', marker='o', linestyle='dashed')
ax[2, 1].plot([0, sync_peak_tau], [0, sync_peak_coherence], color='k', marker='o', linestyle='dashed')

sns.scatterplot(data=df, x='mean_rec_synapses', y='mean_pearson_r', ax=ax[0, 2], color='g')
sns.scatterplot(data=df, x='mean_rec_synapses', y='k_point_1_over_f', ax=ax[1, 2], color='g')
sns.scatterplot(data=df, x='mean_rec_synapses', y='area_over_line', ax=ax[2, 2], color='g')

sns.scatterplot(data=df_normal, x='mean_rec_synapses', y='mean_pearson_r', ax=ax[0, 2], color='k')
sns.scatterplot(data=df_normal, x='mean_rec_synapses', y='k_point_1_over_f', ax=ax[1, 2], color='k')
sns.scatterplot(data=df_normal, x='mean_rec_synapses', y='area_over_line', ax=ax[2, 2], color='k')

