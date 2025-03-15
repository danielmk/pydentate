"""
Simulate the main basket cell ring network.

It produces an .h5 file with the output spikes and several parameters. If 
plotting == True is also performs some analysis and plotting but it does not
save them. The .h5 file is further analyzed in calculate_oscillation_synchrony.py



@author: danielmk
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
from pydentate import net_basket_cell_ring, neuron_tools, oscillations_analysis, spike_to_x
from pydentate.inputs import homogeneous_poisson_process, sigmoid
import os
import matplotlib.pyplot as plt
from scipy.signal import windows, convolve
from scipy.spatial import distance_matrix
from numpy.fft import fft, fftshift, fftfreq
import tables
import matplotlib
import sys
import networkx

dirname = os.path.dirname(__file__)
results_dir = os.path.join(dirname, 'output', 'figure_1_S5_identical_neurons')
if not os.path.isdir(results_dir):
    os.mkdir(results_dir)

seed = 5927

np.random.seed(seed)

plotting = True
gap_resistance = 6e2
gap_delay = 0
dt = 0.1  # In ms
duration = 2  # In s
warmup = 10000  # In ms
input_rate = 10
n_pvbcs = 120
n_input_syns = 100
rec_weight = 7.6e-3
# rec_weight = 0
# pp_bc_weight: 2e-3, 1.9e-3, 1.8e-3, 1.7e-3, 1.6e-3, 1.5e-3, 1.4e-3, 1.3e-3, 1.2e-3, 1.1e-3, 1e-3, 9e-4, 8e-4, 7e-4, 6e-4, 5e-4, 4e-4, 3e-4,, 2e-4, 1e-4
pp_bc_weight = 0
input_current = 0.3  # in nA
input_current_sigma = 0.05  # in nA
fully_connected_on = False
n_rec_syn = 8

gaps_on = False

# Where to search for nrnmech.dll file. Must be adjusted for your machine.
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Generate Input
n_inputs = 192

temporal_patterns = np.array([homogeneous_poisson_process(0.0, duration, input_rate, refractory_period=0.001) for x in range(n_inputs)], dtype=object) * 1000
spatial_patterns = [np.random.choice(range(n_pvbcs), size=n_input_syns, replace=False) for x in range(n_inputs)]

# Create the recurrent connectivity matrix for chemical connections
x = np.linspace(0, np.pi, n_pvbcs)
y = np.sin(x) * 1000
distance_matrix = np.array([np.roll(y, roll) for roll in range(n_pvbcs)])

chem_amplitude = 0.58
chem_center = 141.0
chem_spread = 80

chem_connection_matrix = np.zeros((n_pvbcs, n_pvbcs))
chem_connection_indices = np.array([np.random.choice(n_pvbcs, size=n_rec_syn, replace=False) for x in range(n_pvbcs)])

for row, idc in enumerate(chem_connection_indices):
    # chem_connection_matrix[row, idc] = 1
    chem_connection_matrix[idc, row] = 1

# np.fill_diagonal(chem_connection_matrix, 0)

n_mean_rec_syn = chem_connection_matrix.sum(axis=1).mean()

if fully_connected_on:
    chem_connection_matrix = np.ones(chem_connection_matrix.shape)
    np.fill_diagonal(chem_connection_matrix, 0)

# Create the gap junection connection matrix
# Create the recurrent connectivity matrix for chemical connections
gap_amplitude = 0.773
gap_center = 140.0
gap_spread = 25
probability_matrix = sigmoid(distance_matrix, gap_amplitude, gap_center, gap_spread)

gap_connection_matrix = np.array([[np.random.choice([0, 1], p=[1-x, x]) for x in out] for out in probability_matrix])
np.fill_diagonal(gap_connection_matrix, 0)

"""Create the network"""
nw = net_basket_cell_ring.BasketCellRing(seed+2, temporal_patterns, spatial_patterns, chem_connection_matrix, gap_connection_matrix, rec_weight=rec_weight, n_bcs=n_pvbcs, gap_resistance=6e2, gap_delay=0.0, gap_junctions=gaps_on, pp_bc_weight=pp_bc_weight)  # 7.6e-3

# sys.exit()

"""Create Current Clamp input"""
current_list = []
for cell in nw.populations[0]:
    # rnd_current = np.random.normal(loc=input_current, scale=0)
    current_list.append(input_current)
    cell._current_clamp_soma(amp=input_current, dur=duration*1000, delay=0)

nw.populations[0].voltage_recording(range(120))

"""Run the simulation"""
neuron_tools.run_neuron_simulator(warmup=warmup, t_stop=duration*1000, dt_sim=dt)

"""Save the output spikes"""
times, units = nw.populations[0].get_times_units()

fname = f'pvring_seed_rec-weight_n-pvbcs_gap-resistance_gap-junctions_input-current_input-current-sigma_n-rec-syn_{seed}_{rec_weight}_{n_pvbcs}_{gap_resistance}_{gaps_on}_{input_current}_{input_current_sigma}_{n_rec_syn}.h5'
output_file_path = os.path.join(results_dir, fname)

output_file = tables.File(output_file_path, mode='a')

output_file.create_array('/', 'times', obj=times)
output_file.create_array('/', 'units', obj=units)

parameters_group = output_file.create_group('/', 'parameters')
output_file.create_array('/parameters', 'dt', obj=dt)
output_file.create_array('/parameters', 'seed', obj=seed)
output_file.create_array('/parameters', 'duration', obj=duration)
output_file.create_array('/parameters', 'input_rate', obj=input_rate)
output_file.create_array('/parameters', 'n_pvbcs', obj=n_pvbcs)
output_file.create_array('/parameters', 'n_input_syns', obj=n_input_syns)
output_file.create_array('/parameters', 'n_inputs', obj=n_inputs)
output_file.create_array('/parameters', 'chem_connection_matrix', obj=chem_connection_matrix)
output_file.create_array('/parameters', 'n_rec_synapses', obj=chem_connection_matrix.sum(axis=1))
output_file.create_array('/parameters', 'gap_connection_matrix', obj=gap_connection_matrix)
output_file.create_array('/parameters', 'gap_resistance', obj=gap_resistance)
output_file.create_array('/parameters', 'gap_delay', obj=gap_delay)
output_file.create_array('/parameters', 'input_current', obj=input_current)
output_file.create_array('/parameters', 'input_current_sigma', obj=input_current_sigma)
output_file.create_array('/parameters', 'pp_bc_weight', obj=pp_bc_weight)
output_file.create_array('/parameters', 'gap_junction', obj=gaps_on)
output_file.create_array('/parameters', 'rec_weight', obj=rec_weight)

"""Finish up by closing file"""
output_file.close()

"""Optional analysis/plotting"""
if plotting:
    def exp_decay(t, tau, V):
        '''Exponential decay helper function'''
        return V * np.e ** (-t / tau)


    """Create the synaptic Kernel"""
    decay_tau = 0.0018
    t = np.arange(0,0.020, 0.0001)
    kernel = exp_decay(t, decay_tau, 1)
    kernel = kernel / kernel.sum()
    
    n_timepoints = int((duration * 1000) / dt)
    
    spike_counts = np.unique(times, return_counts=True)
    
    spike_occurences_indices = np.array(spike_counts[0] / dt, dtype=int)
    
    spike_counts_full = np.zeros(n_timepoints+1)
    
    spike_counts_full[spike_occurences_indices] = spike_counts[1]
    
    # hist, bins = np.histogram(output_spiketimes, bins=1000)
    
    sampling_rate = 1 / (dt / 1000)
    # sampling_rate = 5000
    

    plt.rcParams['svg.fonttype'] = 'none'
    
    font = {'family' : 'Arial',
            'size'   : 22}
    matplotlib.rc('font', **font)

    # gap_graph = networkx.from_numpy_matrix(gap_connection_matrix, create_using=networkx.Graph)
    # networkx.draw_circular(gap_graph)
    
    """Coherence Analysis"""
    curr_binary = spike_to_x.units_times_to_binary(np.array(units, dtype=int), times,
                                                   n_units=n_pvbcs,
                                                   dt=dt,
                                                   total_time=duration * 1000 + 1)
    curr_binary = np.array(curr_binary, dtype=int)
    
    curr_binary_convolved = np.array([np.convolve(x, kernel, 'full') for x in curr_binary])
    pr = np.corrcoef(x=curr_binary_convolved)
    ut_idc = np.triu_indices(pr.shape[0], 1)
    mean_pearson_r = np.nanmean(pr[ut_idc])
    
    synaptic_field = curr_binary_convolved.sum(axis=0)
    
    t = np.arange(0, 2001.0, 0.1)
    t_conv_full = np.arange(0, synaptic_field.shape[0] * 0.1, 0.1)

    fig, ax = plt.subplots(2, 1)
    ax[1].plot(t_conv_full, synaptic_field)
    ax[1].set_xlabel('Time (ms)')
    ax[1].set_ylabel('Amplitude')
    ax[1].set_xlim(0, 2000)
    ax[0].set_xlim(0, 2000)
    
    fig.tight_layout()
    
    ax[0].scatter(times, units, marker='|')
    
    bin_size_ms = [0.1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35, 40]
    coherence_list = []
    for bs in bin_size_ms:
        coherence = np.nanmean(oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(bs/dt)))
        coherence_list.append(coherence)

    plt.figure()
    plt.plot(bin_size_ms, coherence_list, marker='o')
    plt.ylim((0, 1))
    plt.xlabel("Bin Size (ms)")
    plt.ylabel("Average Coherence")
    
    units_counts = np.unique(units, return_counts=True)
    print(f"Average frequency: {(units_counts[1] / duration).mean()}")
    
    gaussian_sigma = 1 / dt
    gaussian_width = gaussian_sigma * 10
    gaussian_kernel = windows.gaussian(gaussian_width, gaussian_sigma)
    gaussian_kernel = gaussian_kernel / gaussian_kernel.sum()
    
    spike_counts_full_convolved = convolve(spike_counts_full, gaussian_kernel)
    
    """Progress Report Plot"""
    plt.figure()
    plt.scatter(times, units, marker='|', color='k')
    plt.xlim((0, 2000))
    plt.xlabel("Time (ms)")
    plt.ylabel("Unit #")
    
    plt.figure()
    plt.plot(bin_size_ms, coherence_list, marker='o', color='k')
    plt.ylim((0, 1))
    plt.xlabel("Bin Size (ms)")
    plt.ylabel("Average Coherence")
    plt.xlim((0, 30))

    """Shuffling Oscillation Analysis"""
    # shuffle = oscillations_analysis.shuffling_coherence(times, dt, duration)
    
    """Linearity Analysis"""
    units_counts = np.unique(units, return_counts=True)
    avg_freq = (units_counts[1] / duration).mean()
    
    avg_freq_tau = ((1/avg_freq) * 1000)  # ms

    """Calculate Coherence Measures"""   
    area_over_line, area_over_line_normalized = oscillations_analysis.linearity_analysis(curr_binary, dt, duration, n_points=30)
    
    def exp_decay(t, tau, V):
        '''Exponential decay helper function'''
        return V * np.e ** (-t / tau)

    delta_t_point_1 = (0.1 / avg_freq) * 1000

    k_point_one_over_f = np.nanmean(oscillations_analysis.pairwise_coherence(curr_binary, bin_size=int(delta_t_point_1 / dt)))

    """Create the synaptic Kernel"""   
    curr_binary = np.array(curr_binary, dtype=int)
    
    curr_binary_convolved = np.array([np.convolve(x, kernel, 'full') for x in curr_binary])
    pr = np.corrcoef(x=curr_binary_convolved)
    ut_idc = np.triu_indices(pr.shape[0], 1)
    mean_pearson_r = np.nanmean(pr[ut_idc])
    
    
    
    
    
    

