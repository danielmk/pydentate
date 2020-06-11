# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""


from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
import net_tunedrev
import os
import argparse

# Handle command line inputs with argparse
parser = argparse.ArgumentParser(description='Spatial inhibition paradigm')
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
parser.add_argument('-n_cells',
                    type=int,
                    help='number of cells to stimulate',
                    default=60,
                    dest='n_cells')
parser.add_argument('-cellstomeasure',
                    nargs=3,
                    type=int,
                    help='start stop range for GCs that are SEClamp measured',
                    default=[0, 2000, 50],
                    dest='cellstomeasure')

args = parser.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
n_cells = args.n_cells
cells_to_measure = np.arange(args.cellstomeasure[0],
                             args.cellstomeasure[1],
                             args.cellstomeasure[2])

# Locate and load the nrnmech.dll file. Must to be adjusted for your machine.
dll_files = [("C:\\Users\\DanielM\\Repos\\models_dentate\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll"),
             "C:\\Users\\daniel\\repos\\nrnmech.dll",
             ("C:\\Users\\Holger\\danielm\\models_dentate\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll"),
             ("C:\\Users\\Daniel\\repos\\"
              "dentate_gyrus_Santhakumar2005_and_Yim_patterns\\"
              "dentategyrusnet2005\\nrnmech.dll")]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

# Setup specs for stimulation
stim_pool = 500  # Size of the pool from which stimulated cells are chosen
stim_loc = int(2000 / 2.0 - stim_pool / 2.0)
stim_amp = 1.5
stim_dur = 10
stim_delay = 50

dt = 0.01
for run in runs:
    # Create a standard networks and add the stimulation
    nw = net_tunedrev.TunedNetwork(seed=10000+run)
    np.random.seed(10000 + run)

    # Make sure we are not stimulating a cell we measure
    stim_cells = np.random.choice(range(stim_loc, stim_loc+stim_pool),
                                  n_cells, replace=False)
    while np.intersect1d(stim_cells, cells_to_measure).any():
        stim_cells = np.random.choice(range(stim_loc, stim_loc+stim_pool),
                                      n_cells, replace=False)
    print("Done intersecting")

    nw.populations[0].current_clamp_range(stim_cells,
                                          amp=stim_amp,
                                          dur=stim_dur,
                                          delay=stim_delay)

    nw.populations[0].SEClamp(cells_to_measure, dur1=100, rs=1)
    nw.populations[0].voltage_recording(stim_cells)

    """Initialization for -2000 to -100"""
    h.cvode.active(0)
    h.steps_per_ms = 1.0/dt
    h.finitialize(-60)
    h.t = -2000
    h.secondorder = 0
    h.dt = 10
    while h.t < -100:
        h.fadvance()
    h.secondorder = 2
    h.t = 0
    h.dt = 0.01

    """Setup run control for -100 to 1500"""
    h.frecord_init()  # Necessary after changing t to restart the vectors
    while h.t < 100:
        h.fadvance()

    spk_plot_fname = (str(nw) + "_spikeplot_run_n-cells_cellstomeasure_" +
                      str(run) + '_' + str(n_cells) + '_' +
                      str(args.cellstomeasure))
    data_file_name = (str(nw) + "_data_run_n-cells_cellstomeasure_" +
                      str(run) + '_' + str(n_cells) + '_' +
                      str(args.cellstomeasure))
    spike_plot = nw.plot_aps()
    nw.save_ap_fig(spike_plot, directory=savedir,
                   file_name=spk_plot_fname)
    nw.shelve_network(directory=savedir, file_name=data_file_name)
