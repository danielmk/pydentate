# -*- coding: utf-8 -*-
"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run

@author: DanielM
"""

from neuron import h
import numpy as np
import net_tuned
import matplotlib.pyplot as plt

# Office PC
h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")
# Home PC
#h.nrn_load_dll("C:\\Users\\daniel\\repos\\nrnmech.dll")

# Setup specs for stimulation
n_cells = 100  # Number of cells that are stimulated
stim_pool = 150  # Size of the pool from which stimulated cells are chosen
stim_location = int(2000 / 2.0 - stim_pool / 2.0)
stim_amp = 1.5
stim_dur = 15
stim_delay = 50

# Setup specs for measurements
cells_to_measure = np.arange(0, 2000, 50)

save_dir = "C:\\Users\\DanielM\\Repos\\pyDentate\\paradigm_spatial_inhibition_saves_2018-03-23"

for run in range(2,10):
    # Create a standard networks and add the stimulation
    nw = net_tuned.TunedNetwork(seed=10000+run)
    np.random.seed(10000 + run)
    
    # Make sure we are not stimulating a cell we measure
    stim_cells = np.random.choice(range(stim_location, stim_location+stim_pool), n_cells, replace = False)
    while np.intersect1d(stim_cells, cells_to_measure).any():
        stim_cells = np.random.choice(range(stim_location, stim_location+stim_pool), n_cells, replace = False)
    print("Done intersecting")
    
    nw.populations[0].current_clamp_range(stim_cells,
                                          amp=stim_amp,
                                          dur=stim_dur,
                                          delay=stim_delay)

    nw.populations[0].SEClamp(cells_to_measure, dur1 = 100, rs=1)
    nw.populations[0].voltage_recording(stim_cells)

    """Initialization for -2000 to -100"""
    #print("Running trial " + str(trial))
    h.cvode.active(0)
    dt = 0.01
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
    while h.t < 100:
        h.fadvance()

    spike_plot = nw.plot_aps()
    spike_plot_file_name = "run_" + str(run) + "_spike_plot"
    nw.save_ap_fig(spike_plot, directory = save_dir, file_name = spike_plot_file_name)
    data_file_name = "run_" + str(run) + "_data"
    nw.shelve_network(directory = save_dir, file_name = data_file_name)
    
    # Calculate spatial IPSC plot
    sampling_period = h.dt
    bl_times = np.array([40, 50])  # in ms
    IPSC_times = np.array([50, 90])  # in ms
    bl_dtps = bl_times / sampling_period
    IPSC_dtps = IPSC_times / sampling_period
    
    IPSCs = []
    for cell_i in nw.populations[0].VClamps_i:
        trace = cell_i.as_numpy()
        bl = trace[int(bl_dtps[0]):int(bl_dtps[1])].mean()
        peak_IPSC = trace[int(IPSC_dtps[0]):int(IPSC_dtps[1])].max()
        IPSCs.append(peak_IPSC - bl)
    spatial_plot = plt.figure()
    plt.plot(cells_to_measure, IPSCs)
    plt.xlabel("Granule Cell #")
    plt.ylabel("Peak IPSC (nA)")
    full_file_path = save_dir + '\\' + 'run_' + str(run) + '_spatial_IPSC_plot'
    spatial_plot.savefig(full_file_path + ".pdf", dpi = 300, format ='pdf')
    spatial_plot.savefig(full_file_path + ".eps", dpi = 300, format ='eps')

    

