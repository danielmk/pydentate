"""
Simulate the main basket cell ring network.

It produces an .h5 file with the output spikes and several parameters. If 
plotting == True is also performs some analysis and plotting but it does not
save them. The .h5 file is further analyzed in calculate_oscillation_synchrony.py

@author: danielmk
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from neuron import h, gui  # gui necessary for some parameters to h namespace
from scipy.signal import windows, convolve
from scipy.spatial import distance_matrix
from numpy.fft import fft, fftshift, fftfreq
from pydentate import net_basket_cell_ring, neuron_tools, oscillations_analysis, spike_to_x
from pydentate.inputs import homogeneous_poisson_process, sigmoid

# Directories
dirname = os.path.dirname(__file__)
results_dir = os.path.join(dirname, 'output')
if not os.path.isdir(results_dir):
    os.mkdir(results_dir)

# Constants
SEED = 467
GAP_RESISTANCE = 600  # in Ohms
GAP_DELAY = 0  # in ms
DT = 0.1  # in ms
DURATION = 2  # in s
WARMUP = 2000  # in ms
INPUT_RATE = 10
N_PVBCS = 120
N_INPUT_SYNS = 100
REC_WEIGHT = 7.6e-3
PP_BC_WEIGHT = 0
INPUT_CURRENT = 0.3  # in nA
INPUT_CURRENT_SIGMA = 0.05  # in nA
FULLY_CONNECTED_ON = False
GAPS_ON = True

# Seed for reproducibility
np.random.seed(SEED)

# Load compiled mechanisms
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Generate Input
N_INPUTS = 192

temporal_patterns = np.array([homogeneous_poisson_process(0.0, DURATION, INPUT_RATE, refractory_period=0.001) for _ in range(N_INPUTS)], dtype=object) * 1000
spatial_patterns = [np.random.choice(range(N_PVBCS), size=N_INPUT_SYNS, replace=False) for _ in range(N_INPUTS)]

# Create the recurrent connectivity matrix for chemical connections
x = np.linspace(0, np.pi, N_PVBCS)
y = np.sin(x) * 1000
dist_matrix = np.array([np.roll(y, roll) for roll in range(N_PVBCS)])

chem_amplitude = 0.58
chem_center = 141.0
chem_spread = 80

probability_matrix = sigmoid(dist_matrix, chem_amplitude, chem_center, chem_spread)

chem_connection_matrix = np.array([[np.random.choice([0, 1], p=[1-x, x]) for x in out] for out in probability_matrix])
np.fill_diagonal(chem_connection_matrix, 0)

n_mean_rec_syn = chem_connection_matrix.sum(axis=1).mean()

if FULLY_CONNECTED_ON:
    chem_connection_matrix = np.ones(chem_connection_matrix.shape)
    np.fill_diagonal(chem_connection_matrix, 0)

# Create the gap junction connection matrix
gap_amplitude = 0.773
gap_center = 140.0
gap_spread = 25

probability_matrix = sigmoid(dist_matrix, gap_amplitude, gap_center, gap_spread)

gap_connection_matrix = np.array([[np.random.choice([0, 1], p=[1-x, x]) for x in out] for out in probability_matrix])
np.fill_diagonal(gap_connection_matrix, 0)

# Graphs for chemical and gap junction connections
chem_graph = nx.DiGraph(chem_connection_matrix)
gap_graph = nx.DiGraph(gap_connection_matrix)

colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']

# Plotting chemical connections
nx.draw_networkx_nodes(chem_graph, pos=nx.circular_layout(chem_graph), node_size=10, node_color='k', alpha=1.0)
nx.draw_networkx_edges(chem_graph, pos=nx.circular_layout(chem_graph), edge_color='#8A963F', alpha=0.4)

# Plotting gap junction connections
nx.draw_networkx_nodes(gap_graph, pos=nx.circular_layout(gap_graph), node_size=10, node_color='k', alpha=1.0)
nx.draw_networkx_edges(gap_graph, pos=nx.circular_layout(gap_graph), edge_color='#EA521C', alpha=0.4)

plt.show()
