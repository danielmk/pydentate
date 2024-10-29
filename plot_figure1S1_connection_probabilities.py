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
import matplotlib as mpl


dirname = os.path.dirname(__file__)
results_dir = os.path.join(dirname, 'output')
if not os.path.isdir(results_dir):
    os.mkdir(results_dir)

seed = 467

np.random.seed(seed)

plotting = True
gap_resistance = 6e2
gap_delay = 0
dt = 0.1  # In ms
duration = 2  # In s
warmup = 2000  # In ms
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

gaps_on = True

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

probability_matrix = sigmoid(distance_matrix, chem_amplitude, chem_center, chem_spread)

chem_connection_matrix = np.array([[np.random.choice([0, 1], p=[1-x, x]) for x in out] for out in probability_matrix])
np.fill_diagonal(chem_connection_matrix, 0)

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

colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']
plt.rcParams['svg.fonttype'] = 'none'
font = {'family' : 'Arial',
        'size'   : 32}
mpl.rc('font', **font)

x = np.arange(0.0, 500.0, 0.1)

y_chem = sigmoid(x, chem_amplitude, chem_center, chem_spread)

y_gap = sigmoid(x, gap_amplitude, gap_center, gap_spread)

fig, ax = plt.subplots(1, 2)

ax[0].plot(x, y_chem, color='#8A963F')
ax[1].plot(x, y_gap, color='#EA521C')

ax[0].set_xlabel("Distance ($\mu$m)")
ax[1].set_xlabel("Distance ($\mu$m)")

ax[0].set_ylabel("Chemical Synapse Connection Probability")
ax[1].set_ylabel("Gap Junction Connection Probability")



for a in ax:
    a.set_xlim((0, 500))
    a.set_ylim((0, 1))


