#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 15:08:21 2019

@author: barisckuru
"""


import numpy as np
import matplotlib.pyplot as plt
from tmgexp2_simulator import simulate
from peak_finder_tmgexp2 import peakfinder


'''This script simulates the two types of synapse: from granule cell to mossy
cell and from granule cell to interneuron. It also takes normalized data and
convolves it to get rid of stimulation artifacts. Then it plots simulated G
response with initial hand-tuned taus, optimized taus, and plots peaks of
data for diff freqs'''

gcin_data_path = "/home/can/Downloads/gc_to_in/"
gcin_peaks = peakfinder(gcin_data_path)[1]  # [1] to take normalized values
gcmc_data_path = "/home/can/Downloads/gc_to_mc/"
gcmc_peaks = peakfinder(gcmc_data_path)[1]

x3 = np.array([4.11189527e+03, 8.78187450e-03, 8.05591091e-02])  # gc to in opt
x1 = np.array([3.82825788e+03, 3.23943913e+01, 7.80077153e-02])  # mc to gc opt
x0 = np.array([500, 0, 0.1])  # initial hand-tuned

"CONVOLUTION"
# low pass filter, convolution with 30 datapoints
filt = np.ones(30)/30
gcmc_1Hz = np.convolve(gcmc_peaks[0], filt)
gcmc_10Hz = np.convolve(gcmc_peaks[1], filt)
gcmc_30Hz = np.convolve(gcmc_peaks[2], filt)
gcmc_50Hz = np.convolve(gcmc_peaks[3], filt)

gcin_1Hz = np.convolve(gcin_peaks[0], filt)
gcin_10Hz = np.convolve(gcin_peaks[1], filt)
gcin_30Hz = np.convolve(gcin_peaks[2], filt)
gcin_50Hz = np.convolve(gcin_peaks[3], filt)

'Plotting'


# set_xlim was used to cut the interval we want
# HEX color codes were obtained from colorbrewer2.org website
# codes selected are quantiative and colorblindsafe


fig, axes = plt.subplots(2, 4, sharey=True)

pre_1 = simulate(x0[0], x0[1], 1, x0[2], 0.04975)[1]
opt_1 = simulate(x1[0], x1[1], 1, x1[2], 0.04975)[1]
axes[0, 0].plot(pre_1, color='#1b9e77', linewidth=1)
axes[0, 0].plot(opt_1, color='#d95f02', linewidth=1)
axes[0, 0].plot(gcmc_1Hz, color='#7570b3', linewidth=1)
axes[0, 0].legend(("Initial", "Optimized ", "Data 1Hz"))
axes[0, 0].set_xlim(0, 200000)

pre_10 = simulate(x0[0], x0[1], 10, x0[2], 0.05)[1]
opt_10 = simulate(x1[0], x1[1], 10, x1[2], 0.05)[1]
axes[0, 1].plot(pre_10, color='#1b9e77', linewidth=1)
axes[0, 1].plot(opt_10, color='#d95f02', linewidth=1)
axes[0, 1].plot(gcmc_10Hz, color='#7570b3', linewidth=1)
axes[0, 1].legend(("Initial", "Optimized ", "Data 10Hz"))
axes[0, 1].set_xlim(9000, 30000)

pre_30 = simulate(x0[0], x0[1], 30, x0[2], 0.05)[1]
opt_30 = simulate(x1[0], x1[1], 30, x1[2], 0.05)[1]
axes[0, 2].plot(pre_30, color='#1b9e77', linewidth=1)
axes[0, 2].plot(opt_30, color='#d95f02', linewidth=1)
axes[0, 2].plot(gcmc_30Hz, color='#7570b3', linewidth=1)
axes[0, 2].legend(("Initial", "Optimized ", "Data 30Hz"))
axes[0, 2].set_xlim(9000, 19000)

pre_50 = simulate(x0[0], x0[1], 50, x0[2], 0.05)[1]
opt_50 = simulate(x1[0], x1[1], 50, x1[2], 0.05)[1]
axes[0, 3].plot(pre_50, color='#1b9e77', linewidth=1)
axes[0, 3].plot(opt_50, color='#d95f02', linewidth=1)
axes[0, 3].plot(gcmc_50Hz, color='#7570b3', linewidth=1)
axes[0, 3].legend(("Initial", "Optimized", "Data 50Hz"))
axes[0, 3].set_xlim(9000, 15000)

axes[0, 1].set_title("Granule cell to Mossy cell")


gcin_pre_1 = simulate(x0[0], x0[1], 1, x0[2], 0.04975)[1]
gcin_opt_1 = simulate(x3[0], x3[1], 1, x3[2], 0.04975)[1]
axes[1, 0].plot(gcin_pre_1, color='#1b9e77', linewidth=1)
axes[1, 0].plot(gcin_opt_1, color='#d95f02', linewidth=1)
axes[1, 0].plot(gcin_1Hz, color='#7570b3', linewidth=1)
axes[1, 0].legend(("Initial", "Optimized ", "Data 1Hz"))
axes[1, 0].set_xlim(0, 200000)

gcin_pre_10 = simulate(x0[0], x0[1], 10, x0[2], 0.05)[1]
gcin_opt_10 = simulate(x3[0], x3[1], 10, x3[2], 0.05)[1]
axes[1, 1].plot(gcin_pre_10, color='#1b9e77', linewidth=1)
axes[1, 1].plot(gcin_opt_10, color='#d95f02', linewidth=1)
axes[1, 1].plot(gcin_10Hz, color='#7570b3', linewidth=1)
axes[1, 1].legend(("Initial", "Optimized ", "Data 10Hz"))
axes[1, 1].set_xlim(9000, 30000)

gcin_pre_30 = simulate(x0[0], x0[1], 30, x0[2], 0.05)[1]
gcin_opt_30 = simulate(x3[0], x3[1], 30, x3[2], 0.05)[1]
axes[1, 2].plot(gcin_pre_30, color='#1b9e77', linewidth=1)
axes[1, 2].plot(gcin_opt_30, color='#d95f02', linewidth=1)
axes[1, 2].plot(gcin_30Hz, color='#7570b3', linewidth=1)
axes[1, 2].legend(("Initial", "Optimized ", "Data 30Hz"))
axes[1, 2].set_xlim(9000, 19000)

gcin_pre_50 = simulate(x0[0], x0[1], 50, x0[2], 0.05)[1]
gcin_opt_50 = simulate(x3[0], x3[1], 50, x3[2], 0.05)[1]
axes[1, 3].plot(gcin_pre_50, color='#1b9e77', linewidth=1)
axes[1, 3].plot(gcin_opt_50, color='#d95f02', linewidth=1)
axes[1, 3].plot(gcin_50Hz, color='#7570b3', linewidth=1)
axes[1, 3].legend(("Initial", "Optimized", "Data 50Hz"))
axes[1, 3].set_xlim(9000, 15000)

axes[1, 1].set_title("Granule cell to Interneuron")
axes[1, 0].set_xlabel("Data points")
axes[0, 0].set_ylabel("Normalized G")
axes[1, 0].set_ylabel("Normalized G")
