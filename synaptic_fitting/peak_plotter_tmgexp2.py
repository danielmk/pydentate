#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:08:21 2019

@author: can
"""

import numpy as np
import matplotlib.pyplot as plt
from tmgexp2_simulator import simulate
from peak_finder_tmgexp2 import peakfinder

'''This script simulates the two types of synapse: from granule cell to mossy
cell and from granule cell to interneuron and takes the peaks of signals with
sems. Then it plots simulation peaks with initial hand-tuned taus, optimized
taus and plots peaks of data for diff freqs'''


gcin_data_path = "/home/can/Downloads/gc_to_in/"
gcin_peaks = peakfinder(gcin_data_path)[0]  # [1] to take normalized values
# number of cells used in average data, required to calc SEM out of SD
gcin_sem = peakfinder(gcin_data_path)[2]/np.sqrt(9)  # 9 cells
gcmc_data_path = "/home/can/Downloads/gc_to_mc/"
gcmc_peaks = peakfinder(gcmc_data_path)[0]
gcmc_sem = peakfinder(gcmc_data_path)[2]/np.sqrt(10)  # 10 cells

gcmc_1Hz = gcmc_peaks[0]
gcmc_10Hz = gcmc_peaks[1]
gcmc_30Hz = gcmc_peaks[2]
gcmc_50Hz = gcmc_peaks[3]
gcmc_sem_1Hz = gcmc_sem[0]
gcmc_sem_10Hz = gcmc_sem[1]
gcmc_sem_30Hz = gcmc_sem[2]
gcmc_sem_50Hz = gcmc_sem[3]


gcin_1Hz = gcin_peaks[0]
gcin_10Hz = gcin_peaks[1]
gcin_30Hz = gcin_peaks[2]
gcin_50Hz = gcin_peaks[3]
gcin_sem_1Hz = gcin_sem[0]
gcin_sem_10Hz = gcin_sem[1]
gcin_sem_30Hz = gcin_sem[2]
gcin_sem_50Hz = gcin_sem[3]


# Optimized values for tau facilitation, tau recovery and u0 respectively
x1 = np.array([4.11189527e+03, 8.78187450e-03, 8.05591091e-02])  # gc to in opt
x3 = np.array([3.82825788e+03, 3.23943913e+01, 7.80077153e-02])  # gc to mc opt
x0 = np.array([500, 0, 0.1])  # initial hand-tuned


'Plotting'

plt.rcParams.update({'errorbar.capsize': 2})
fig, axes = plt.subplots(2, 4, sharey=True)
a = np.arange(1, 11)

pre_1 = simulate(x0[0], x0[1], 1, x0[2])[0]
opt_1 = simulate(x1[0], x1[1], 1, x1[2])[0]
axes[0, 0].plot(a, pre_1, color='#1b9e77', marker='o', markersize=3)
axes[0, 0].plot(a, opt_1, color='#d95f02', marker='o', markersize=3)
axes[0, 0].errorbar(a, gcmc_1Hz, gcmc_sem_1Hz, color='#7570b3', marker='o', markersize=3)
axes[0, 0].legend(("Initial", "Optimized ", "Data 1Hz"))

pre_10 = simulate(x0[0], x0[1], 10, x0[2])[0]
opt_10 = simulate(x1[0], x1[1], 10, x1[2])[0]
axes[0, 1].plot(a, pre_10, color='#1b9e77', marker='o', markersize=3)
axes[0, 1].plot(a, opt_10, color='#d95f02', marker='o', markersize=3)
axes[0, 1].errorbar(a, gcmc_10Hz, gcmc_sem_10Hz, color='#7570b3', marker='o', markersize=3)
axes[0, 1].legend(("Initial", "Optimized ", "Data 10Hz"))

pre_30 = simulate(x0[0], x0[1], 30, x0[2])[0]
opt_30 = simulate(x1[0], x1[1], 30, x1[2])[0]
axes[0, 2].plot(a, pre_30, color='#1b9e77', marker='o', markersize=3)
axes[0, 2].plot(a, opt_30, color='#d95f02', marker='o', markersize=3)
axes[0, 2].errorbar(a, gcmc_30Hz, gcmc_sem_30Hz, color='#7570b3', marker='o', markersize=3)
axes[0, 2].legend(("Initial", "Optimized ", "Data 30Hz"))

pre_50 = simulate(x0[0], x0[1], 50, x0[2])[0]
opt_50 = simulate(x1[0], x1[1], 50, x1[2])[0]
axes[0, 3].plot(a, pre_50, color='#1b9e77', marker='o', markersize=3)
axes[0, 3].plot(a, opt_50, color='#d95f02', marker='o', markersize=3)
axes[0, 3].errorbar(a, gcmc_50Hz, gcmc_sem_50Hz, color='#7570b3', marker='o', markersize=3)
axes[0, 3].legend(("Initial", "Optimized", "Data 50Hz"))

axes[0,1].set_title("Granule cell to Mossy cell")


gcin_pre_1 = simulate(x0[0], x0[1], 1, x0[2])[0]
gcin_opt_1 = simulate(x3[0], x3[1], 1, x3[2])[0]
axes[1, 0].plot(a, gcin_pre_1, color='#1b9e77', marker='o', markersize=3)
axes[1, 0].plot(a, gcin_opt_1, color='#d95f02', marker='o', markersize=3)
axes[1, 0].errorbar(a, gcin_1Hz, gcin_sem_1Hz, color='#7570b3', marker='o', markersize=3)
axes[1, 0].legend(("Initial", "Optimized ", "Data 1Hz"))

gcin_pre_10 = simulate(x0[0], x0[1], 10, x0[2])[0]
gcin_opt_10 = simulate(x3[0], x3[1], 10, x3[2])[0]
axes[1, 1].plot(a, gcin_pre_10, color='#1b9e77', marker='o', markersize=3)
axes[1, 1].plot(a, gcin_opt_10, color='#d95f02', marker='o', markersize=3)
axes[1, 1].errorbar(a, gcin_10Hz, gcin_sem_10Hz, color='#7570b3', marker='o', markersize=3)
axes[1, 1].legend(("Initial", "Optimized ", "Data 10Hz"))

gcin_pre_30 = simulate(x0[0], x0[1], 30, x0[2])[0]
gcin_opt_30 = simulate(x3[0], x3[1], 30, x3[2])[0]
axes[1, 2].plot(a, gcin_pre_30, color='#1b9e77', marker='o', markersize=3)
axes[1, 2].plot(a, gcin_opt_30, color='#d95f02', marker='o', markersize=3)
axes[1, 2].errorbar(a, gcin_30Hz, gcin_sem_30Hz, color='#7570b3', marker='o', markersize=3)
axes[1, 2].legend(("Initial", "Optimized ", "Data 30Hz"))

gcin_pre_50 = simulate(x0[0], x0[1], 50, x0[2])[0]
gcin_opt_50 = simulate(x3[0], x3[1], 50, x3[2])[0]
axes[1, 3].plot(a, gcin_pre_50, color='#1b9e77', marker='o', markersize=3)
axes[1, 3].plot(a, gcin_opt_50, color='#d95f02', marker='o', markersize=3)
axes[1, 3].errorbar(a, gcin_50Hz, gcin_sem_50Hz, color='#7570b3', marker='o', markersize=3)
axes[1, 3].legend(("Initial", "Optimized", "Data 50Hz"))

axes[1, 1].set_title("Granule cell to Interneuron")
axes[1, 0].set_xlabel("Peaks")
axes[0, 0].set_ylabel("Normalized G")
axes[1, 0].set_ylabel("Normalized G")
