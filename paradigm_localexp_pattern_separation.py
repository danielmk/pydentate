# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h
import numpy as np
import net_tunedrevexp
from burst_generator_inhomogeneous_poisson import inhom_poiss
import os
import argparse
import scipy.stats as stats
import math

# Setup some helper functions
def pos(rad):
    """
    (x,y) position of a point on a circle with axis origin at (0,0)
    and radius 1.
    x = cx + r * cos(rad) -> x = cos(rad)
    y = cy + r * sin(rad) -> y = sin(rad)
    
    Returns a list of tuples that give the point of each radian passed.
    """
    x_arr = list(np.cos(rad))
    y_arr = list(np.sin(rad))
    
    return [(x_arr[idx],y_arr[idx]) for idx in range(len(x_arr))]

def euclidian_dist(p1,p2):
    """ p1 and p2 must both be of len 2 where p1 = (x1,y1); p2 = (x2,y2)"""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
    
    
    

# Handle command line inputs with argparse
parser = argparse.ArgumentParser(description='Local pattern separation paradigm')
parser.add_argument('-runs',
                    nargs=3,
                    type=int,
                    help='start stop range for the range of runs',
                    default=[0, 1, 1],
                    dest='runs')
parser.add_argument('-savedir',
                    type=str,
                    help='complete directory where data is saved',
                    default=os.getcwd(),
                    dest='savedir')
parser.add_argument('-scale',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=500,
                    dest='input_scale')
parser.add_argument('-seed',
                    type=int,
                    help='standard deviation of gaussian distribution',
                    default=0,
                    dest='seed')

args = parser.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
input_scale = args.input_scale
seed = args.seed

# Locate a nrnmech.dll file that has the mechanisms required by the network
# On your own machine you have to add the path to your own file to the list dll_files
dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
            "C:\\Users\\daniel\\Repos\\nrnmech.dll",
            "C:\\Users\\Holger\\danielm\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
            "C:\\Users\\Daniel\\repos\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll"]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

np.random.seed(seed)
# Generate a gaussian probability density function
gauss_gc = stats.expon(loc=0, scale=input_scale)
gauss_bc = stats.expon(loc=0, scale=(input_scale/2000.0)*24)
pdf_gc = gauss_gc.pdf(np.arange(2000))
pdf_gc = pdf_gc/pdf_gc.sum()
pdf_bc = gauss_bc.pdf(np.arange(24))
pdf_bc = pdf_bc/pdf_bc.sum()
# We hold the pdf constant. To randomize the centroid we reslice the GC indices
GC_indices = np.arange(2000, dtype=float)
GC_pop_rads = (GC_indices / 2000.0) * (2*np.pi)
GC_pop_pos = pos(GC_pop_rads)
PP_rads = (np.random.randint(0, 2000, size=400) / 2000.0) * (2*np.pi)
PP_pos =  pos(PP_rads)

PP_to_GCs = []
for curr_PP_pos in PP_pos:
    curr_dist = []
    for post_cell_pos in GC_pop_pos:
        curr_dist.append(euclidian_dist(curr_PP_pos, post_cell_pos))
    sort_idc = np.argsort(curr_dist)
    picked_cells = np.random.choice(sort_idc, 800, replace=True, p = pdf_gc)
    PP_to_GCs.append(picked_cells)
PP_to_GCs = np.array(PP_to_GCs)
# Generate the PP -> BC mapping as above
BC_indices = np.arange(24)
BC_pop_rads = (BC_indices / 24.0) * (2*np.pi)
BC_pop_pos = pos(BC_pop_rads)

PP_to_BCs = []
for curr_PP_pos in PP_pos:
    curr_dist = []
    for post_cell_pos in BC_pop_pos:
        curr_dist.append(euclidian_dist(curr_PP_pos, post_cell_pos))
    sort_idc = np.argsort(curr_dist)
    picked_cells = np.random.choice(sort_idc, 10, replace=True, p = pdf_bc)
    PP_to_BCs.append(picked_cells)

PP_to_BCs = np.array(PP_to_BCs)

# Generate temporal patterns for the 100 PP inputs
np.random.seed(seed)
temporal_patterns = inhom_poiss()

# Start the runs of the model
for run in runs:
    nw = net_tunedrevexp.TunedNetwork(seed, temporal_patterns[0+run:24+run],
                                PP_to_GCs[0+run:24+run],
                                PP_to_BCs[0+run:24+run])

    # Attach voltage recordings to all cells
    nw.populations[0].voltage_recording(range(2000))
    nw.populations[1].voltage_recording(range(60))
    nw.populations[2].voltage_recording(range(24))
    nw.populations[3].voltage_recording(range(24))
    # Run the model
    """Initialization for -2000 to -100"""
    h.cvode.active(0)
    dt = 0.1
    h.steps_per_ms = 1.0/dt
    h.tstop = 1500
    h.finitialize(-60)
    h.t = -2000
    h.secondorder = 0
    h.dt = 10
    while h.t < -100:
        h.fadvance()
        print(h.t)

    h.secondorder = 2
    h.t = 0
    h.dt = 0.1

    """Setup run control for -100 to 1500"""
    h.frecord_init() # Necessary after changing t to restart the vectors
    
    while h.t < 600:
        h.fadvance()
    print("Done Running")

    tuned_save_file_name = str(nw) + '_data_paradigm_local-pattern-separation_run_scale_seed_' + str(run).zfill(3) + '_' + str(input_scale).zfill(3) + '_' + str(seed)
    nw.shelve_network(savedir, tuned_save_file_name)

    fig = nw.plot_aps(time=600)
    tuned_fig_file_name = str(nw) + '_spike-plot_paradigm_local-pattern-separation_run_scale_seed_' + str(run).zfill(3) + '_' + str(input_scale).zfill(3) + '_' + str(seed)
    nw.save_ap_fig(fig, savedir, tuned_fig_file_name)
