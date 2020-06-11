# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""


from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
import net_nonfacilitatingrev
import os
import argparse

# Handle command line inputs
parser = argparse.ArgumentParser(description='Frequency inhibition paradigm')
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
parser.add_argument('-interval',
                    type=int,
                    help='stimulus interval in ms',
                    default=1000,
                    dest='interval')
parser.add_argument('-n_cells',
                    type=int,
                    help='number of cells to stimulate',
                    default=40,
                    dest='n_cells')

args = parser.parse_args()
runs = range(args.runs[0], args.runs[1], args.runs[2])
savedir = args.savedir
interval = args.interval
n_cells = args.n_cells

# Where to search for nrnmech.dll file. Must be adjusted for your machine.
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
stim_pool = 150  # Size of the pool from which stimulated cells are chosen
stim_location = int(2000 / 2.0 - stim_pool / 2.0)  # Global because 2000
stim_amp = 1
stim_dur = 5
stim_delay = 100

# Setup specs for measurements
cells_to_measure = np.arange(0, 2000, 100)

dt = 0.01
for run in runs:
    # Create a networks
    nw = net_nonfacilitatingrev.TunedNetwork(seed=10000+run)
    np.random.seed(10000 + run)

    # Make sure we are not stimulating a cell we measure
    stim_cells = np.random.choice(2000, n_cells)
    gcs_to_measure = np.random.choice(2000, 20)
    while np.intersect1d(stim_cells, gcs_to_measure).any():
        gcs_to_measure = np.random.choice(2000, 20)

    for x in range(10):
        nw.populations[0].current_clamp_range(stim_cells,
                                              amp=stim_amp,
                                              dur=stim_dur,
                                              delay=100+interval*x)

    nw.populations[0].SEClamp(gcs_to_measure, dur1=100+200+interval*10, rs=1)
    nw.populations[0].voltage_recording(range(0, 2000, 100))
    nw.populations[1].voltage_recording(range(0, 60, 2))
    nw.populations[2].voltage_recording(range(24))
    nw.populations[3].voltage_recording(range(24))
    mc_to_measure = np.random.choice(60, 1)
    bc_to_measure = np.random.choice(24, 1)
    hc_to_measure = np.random.choice(24, 1)
    nw.populations[1].SEClamp(mc_to_measure, dur1=100+200+interval*10,
                              rs=1, amp1=-70)
    nw.populations[2].SEClamp(bc_to_measure, dur1=100+200+interval*10,
                              rs=1, amp1=-70)
    nw.populations[3].SEClamp(hc_to_measure, dur1=100+200+interval*10,
                              rs=1, amp1=-70)

    """Initialization for -2000 to -100"""
    h.cvode.active(0)
    h.steps_per_ms = 1.0/dt
    h.tstop = 1500
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
    while h.t < 100+200+interval*10:
        h.fadvance()

    spike_plot = nw.plot_aps(time=100+200+interval*10)
    spike_plot_file_name = (str(nw) + "_spikeplot_run_interval_n-cells_" +
                            str(run) + '_' + str(interval) + '_' +
                            str(n_cells))
    nw.save_ap_fig(spike_plot, directory=savedir,
                   file_name=spike_plot_file_name)
    data_file_name = (str(nw) + "_data_run_interval_n-cells_" + str(run) +
                      '_' + str(interval) + '_' + str(n_cells))
    nw.shelve_network(directory=savedir, file_name=data_file_name)
