# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 11:54:22 2019

@author: Daniel
"""

from tmgexp2syn_fitter import simulate, peakfinder
import numpy as np
import matplotlib.pyplot as plt
from mossycell import MossyCell
from basketcell import BasketCell

freqs = [1,10,30,50]
gc_mc_sim_initial = []
gc_mc_sim_optimized = []
gc_in_sim_initial = []
gc_in_sim_optimized = []

for stim_freq in freqs:
    curr_gc_mc_initial = simulate(  3000,
                                    25,
                                    stim_freq,
                                    MossyCell,
                                    tau_1=0.5,
                                    tau_2=6.2,
                                    u0=0.1,
                                    sampling_interval=0.5,
                                    v_clamp=-70.42,
                                    sec='proxd')
    
    curr_gc_mc_optimized = simulate(3.73987038e+03,
                                    3.96206731e+01,
                                    stim_freq,
                                    MossyCell,
                                    tau_1=0.5,
                                    tau_2=6.2,
                                    u0=5.94309891e-02,
                                    sampling_interval=0.5,
                                    v_clamp=-70.42,
                                    sec='proxd')
    curr_gc_in_initial = simulate(  3000,
                                    10,
                                    stim_freq,
                                    MossyCell,
                                    tau_1=0.5,
                                    tau_2=6.2,
                                    u0=0.1,
                                    sampling_interval=0.5,
                                    v_clamp=-70.42,
                                    sec='proxd')
    
    curr_gc_in_optimized = simulate(3.52671493e+03,
                                    1.33249553e+01,
                                    stim_freq,
                                    MossyCell,
                                    tau_1=0.5,
                                    tau_2=6.2,
                                    u0=5.77352002e-02,
                                    sampling_interval=0.5,
                                    v_clamp=-70.42,
                                    sec='proxd')
    gc_mc_sim_initial.append(curr_gc_mc_initial[0])
    gc_mc_sim_optimized.append(curr_gc_mc_optimized[0])
    gc_in_sim_initial.append(curr_gc_in_initial[0])
    gc_in_sim_optimized.append(curr_gc_in_optimized[0])

data_path = (
    "C:\\Users\\Daniel\\Dropbox\\02_MEC Project\\003_Antidromic "
    "electrically evoked EPSCs in hilar cells, voltage clamp\\"
)

gc_to_in = np.load(data_path + "gc_to_in_data_full.npz")
gc_to_mc = np.load(data_path + "gc_to_mc_data_full.npz")

x = range(1,11)
fig, axes = plt.subplots(2, 4)
csiminit = '#1b9e77'
csimopt = '#d95f02'
cdata = '#7570b3'
legend = ("Initial", "Optimized ", "Data")
axes[0,0].errorbar(x, gc_to_mc['peaks_norm_01'], yerr=gc_to_mc['peaks_sem_01'], marker='o', color = cdata,capsize=2, capthick=1)
axes[0,0].plot(x, gc_mc_sim_initial[0], marker='o', color = csiminit)
axes[0,0].plot(x, gc_mc_sim_optimized[0], marker='o', color = csimopt)
axes[0,0].legend(legend)

axes[0,1].errorbar(x, gc_to_mc['peaks_norm_10'], yerr=gc_to_mc['peaks_sem_10'], marker='o', color = cdata,capsize=2, capthick=1)
axes[0,1].plot(x, gc_mc_sim_initial[1], marker='o', color = csiminit)
axes[0,1].plot(x, gc_mc_sim_optimized[1], marker='o', color = csimopt)

axes[0,2].errorbar(x, gc_to_mc['peaks_norm_30'], yerr=gc_to_mc['peaks_sem_30'], marker='o', color = cdata,capsize=2, capthick=1)
axes[0,2].plot(x, gc_mc_sim_initial[2], marker='o', color = csiminit)
axes[0,2].plot(x, gc_mc_sim_optimized[2], marker='o', color = csimopt)

axes[0,3].errorbar(x, gc_to_mc['peaks_norm_50'], yerr=gc_to_mc['peaks_sem_50'], marker='o', color = cdata,capsize=2, capthick=1)
axes[0,3].plot(x, gc_mc_sim_initial[3], marker='o', color = csiminit)
axes[0,3].plot(x, gc_mc_sim_optimized[3], marker='o', color = csimopt)

axes[1,0].errorbar(x, gc_to_in['peaks_norm_01'], yerr=gc_to_in['peaks_sem_01'], marker='o', color = cdata,capsize=2, capthick=1)
axes[1,0].plot(x, gc_in_sim_initial[0], marker='o', color = csiminit)
axes[1,0].plot(x, gc_in_sim_optimized[0], marker='o', color = csimopt)
axes[1,0].legend(legend)

axes[1,1].errorbar(x, gc_to_in['peaks_norm_10'], yerr=gc_to_in['peaks_sem_10'], marker='o', color = cdata,capsize=2, capthick=1)
axes[1,1].plot(x, gc_in_sim_initial[1], marker='o', color = csiminit)
axes[1,1].plot(x, gc_in_sim_optimized[1], marker='o', color = csimopt)

axes[1,2].errorbar(x, gc_to_in['peaks_norm_30'], yerr=gc_to_in['peaks_sem_30'], marker='o', color = cdata,capsize=2, capthick=1)
axes[1,2].plot(x, gc_in_sim_initial[2], marker='o', color = csiminit)
axes[1,2].plot(x, gc_in_sim_optimized[2], marker='o', color = csimopt)

axes[1,3].errorbar(x, gc_to_in['peaks_norm_50'], yerr=gc_to_in['peaks_sem_50'], marker='o', color = cdata,capsize=2, capthick=1)
axes[1,3].plot(x, gc_in_sim_initial[3], marker='o', color = csiminit)
axes[1,3].plot(x, gc_in_sim_optimized[3], marker='o', color = csimopt)

ylims = ((0,4.5), (0,8))

for idx1, x in enumerate(axes):
    for idx2, y in enumerate(x):
        y.set_ylim(ylims[idx1])
        y.set_xticks(range(1,11))
        y.set_title(str(freqs[idx2])+' Hz')
for x in axes[1]:
    x.set_xlabel("# Stimulus")
    
axes[0,0].set_ylabel("GC to MC\nNormalized Peaks")
axes[1,0].set_ylabel("GC to IN\nNormalized Peaks")

