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


dirname = os.path.dirname(__file__)
data_dir = os.path.join(dirname, 'output', 'slow_gamma_input_high_sync')

all_files = [os.path.join(data_dir, x) for x in os.listdir(data_dir)]

bin_size_ms = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]

neuron_types = {'gc': 10000, 'bc': 120, 'hc': 120, 'mc': 300}

def exp_decay(t, tau, V):
    '''Exponential decay helper function'''
    return V * np.e ** (-t / tau)


"""Create the synaptic Kernel"""
decay_tau = 0.0018
t = np.arange(0,0.020, 0.0001)
kernel = exp_decay(t, decay_tau, 1)

for curr_file in all_files:
    print(curr_file)
    curr_data = tables.open_file(curr_file, mode='a')

    if not hasattr(curr_data.root, 'analysis'):
        curr_data.create_group('/', 'analysis')
    """Do the coherence analysis"""
    if not hasattr(curr_data.root.analysis, 'coherence'):
        curr_data.create_group('/analysis', 'coherence')
        for cell in neuron_types.keys():
            curr_units = eval(f"curr_data.root.units.units_{cell}.read()")
            curr_times = eval(f"curr_data.root.times.times_{cell}.read()")
            curr_dt = curr_data.root.parameters.dt.read()
            curr_units = np.array(curr_units, dtype=int)
            curr_duration = curr_data.root.parameters.duration.read()
    
            curr_binary = spike_to_x.units_times_to_binary(curr_units, curr_times,
                                                           n_units=neuron_types[cell],
                                                           dt=curr_dt,
                                                           total_time=curr_data.root.parameters.duration.read() * 1000 + 1)
            
            curr_binary = np.array(curr_binary, dtype=int)
            
            curr_binary_convolved = np.array([np.convolve(x, kernel, 'full') for x in curr_binary])
            pr = np.corrcoef(x=curr_binary_convolved)
            ut_idc = np.triu_indices(pr.shape[0], 1)
            mean_pearson_r = np.nanmean(pr[ut_idc])
            
            synaptic_field = curr_binary_convolved.sum(axis=0)
            # synaptic_field = np.convolve(curr_binary.sum(axis=0), kernel, 'same')
            
            curr_data.create_array('/analysis/coherence', f'synaptic_field_{cell}', obj=synaptic_field)
            
            curr_binary_convolved = np.array([np.convolve(x, kernel, 'same') for x in curr_binary])
            pr = np.corrcoef(x=curr_binary_convolved)
            ut_idc = np.triu_indices(pr.shape[0], 1)
            mean_pearson_r = np.nanmean(pr[ut_idc])

            mean_coherence_k = []
            for bs in bin_size_ms:
                pairwise_coherence = oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(bs/curr_dt))
                mean_coherence_k.append(np.nanmean(pairwise_coherence))
            curr_data.create_array('/analysis/coherence', f'mean_k_{cell}', obj=mean_coherence_k)
            curr_data.create_array('/analysis/coherence', f'mean_k_dt_{cell}', obj=bin_size_ms)
            curr_data.create_array('/analysis/coherence', f'mean_pearson_r_{cell}', obj=mean_pearson_r)

            """Calculate coherence at frequency dependent delta t"""
            unit_counts = np.unique(curr_units, return_counts=True)
            frequencies = np.array(unit_counts[1] / curr_duration)
            mean_frequency = frequencies.mean()
            delta_t_point_1 = (0.1 / mean_frequency) * 1000
            
            k_point_one_over_f = np.nanmean(oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(delta_t_point_1 / curr_dt)))
            curr_data.create_array('/analysis/coherence', f'delta_t_point_1_{cell}', obj=delta_t_point_1)
            curr_data.create_array('/analysis/coherence', f'k_point_1_over_f_{cell}', obj=k_point_one_over_f)
    
            delta_t_1 = (1 / mean_frequency) * 1000
            k_one_over_f = np.nanmean(oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(delta_t_1 / curr_dt)))
            
            curr_data.create_array('/analysis/coherence', f'delta_t_1_{cell}', obj=delta_t_1)
            curr_data.create_array('/analysis/coherence', f'k_1_over_f_{cell}', obj=k_one_over_f)
            
            k_one_over_f = np.nanmean(oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(delta_t_1 / curr_dt)))
    
            linearity_n_points = 30
            area_over_line, area_over_line_normalized = oscillations_analysis.linearity_analysis(curr_binary, curr_dt, curr_duration, n_points=linearity_n_points)
            
            curr_data.create_array('/analysis/coherence', f'area_over_line_{cell}', obj=area_over_line)
            curr_data.create_array('/analysis/coherence', f'area_over_line_normalized_{cell}', obj=area_over_line_normalized)
            curr_data.create_array('/analysis/coherence', f'linearity_n_points_{cell}', obj=linearity_n_points)

    else:
        Warning(f"In {curr_file} coherence already existed in analysis and was skipped.")

    if not hasattr(curr_data.root.analysis, 'fft'):
        curr_data.create_group('/analysis', 'fft')
        for cell in neuron_types.keys():
            curr_units = eval(f"curr_data.root.units.units_{cell}.read()")
            curr_times = eval(f"curr_data.root.times.times_{cell}.read()")
            curr_dt = curr_data.root.parameters.dt.read()
            curr_duration = curr_data.root.parameters.duration.read()
            
            n_timepoints = int((curr_duration * 1000) / curr_dt)
            
            spike_counts = np.unique(curr_times, return_counts=True)
            
            spike_occurences_indices = np.array(spike_counts[0] / curr_dt, dtype=int)
            
            spike_counts_full = np.zeros(n_timepoints+1)
            
            spike_counts_full[spike_occurences_indices] = spike_counts[1]
            
            gaussian_sigma = 1 / curr_dt
            gaussian_width = gaussian_sigma * 10
            gaussian_kernel = windows.gaussian(gaussian_width, gaussian_sigma)
            gaussian_kernel = gaussian_kernel / gaussian_kernel.sum()
            
            spike_counts_full_convolved = convolve(spike_counts_full, gaussian_kernel)
    
            sampling_rate = 1 / (curr_dt / 1000)
            # sampling_rate = 5000
            
            sp = fftshift(fft(spike_counts_full_convolved))
            freq = fftshift(fftfreq(len(spike_counts_full_convolved), 1/sampling_rate))
            
            amp = np.abs(sp)
            
            curr_data.create_array('/analysis/fft', f'amp_{cell}', obj=amp)
            curr_data.create_array('/analysis/fft', f'freq_{cell}', obj=freq)
            curr_data.create_array('/analysis/fft', f'spike_counts_{cell}', obj=spike_counts_full)
            curr_data.create_array('/analysis/fft', f'spike_counts_convolved_{cell}', obj=spike_counts_full_convolved)

    else:
        Warning(f"In {curr_file} fft already existed in analysis and was skipped.")

    if not hasattr(curr_data.root.analysis, 'frequencies'):
        for cell in neuron_types.keys():
            curr_units = eval(f"curr_data.root.units.units_{cell}.read()")
            curr_times = eval(f"curr_data.root.times.times_{cell}.read()")
            curr_dt = curr_data.root.parameters.dt.read()
            curr_duration = curr_data.root.parameters.duration.read()

            unit_counts = np.unique(curr_units, return_counts=True)
            frequencies = np.array(unit_counts[1] / curr_duration)

            curr_data.create_array('/analysis', f'frequencies_{cell}', obj=frequencies)
    else:
        Warning(f"In {curr_file} frequencies already existed in analysis and was skipped.")

    curr_data.close()

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
