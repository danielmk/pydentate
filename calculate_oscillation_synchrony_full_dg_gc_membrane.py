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
from scipy.signal import windows, convolve, welch, periodogram
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
data_dir = os.path.join(dirname, 'output', 'slow_gamma_input_high_sync_gc_membrane_with_dendrite')

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

    curr_membrane = curr_data.root.gc_membrane.read()
    curr_dendrites = curr_data.root.gc_dendrite.read()

    # sys.exit()

    sampling_rate = 1 / 0.0001
    
    if not hasattr(curr_data.root, 'Pxx_den'):
        f_soma, Pxx_den_soma = periodogram(curr_membrane.mean(axis=0), sampling_rate)
        f_dendrite, Pxx_den_dendrite = periodogram(curr_dendrites.mean(axis=0), sampling_rate)
        # plt.semilogy(f, Pxx_den)
        curr_data.create_array('/', 'frequency', obj=f_soma)
        curr_data.create_array('/', 'Pxx_den', obj=Pxx_den_soma)
        curr_data.create_array('/', 'frequency_dendrites', obj=f_dendrite)
        curr_data.create_array('/', 'Pxx_den_dendrites', obj=Pxx_den_dendrite)
    
    # f_welch, Pxx_den_welch = periodogram(curr_membrane.mean(axis=0), sampling_rate, window=('hamming'))
    
    # curr_data.create_array('/', 'frequency_welch', obj=f_welch)
    # curr_data.create_array('/', 'Pxx_den_welch', obj=Pxx_den_welch)

    if not hasattr(curr_data.root, 'analysis'):
        pass
        # curr_data.create_group('/', 'analysis')
    """Do the coherence analysis"""
    # if not hasattr(curr_data.root.analysis, 'coherence'):
        # curr_data.create_group('/analysis', 'coherence')

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
