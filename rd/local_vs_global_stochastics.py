# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 15:56:08 2018

@author: Daniel
"""
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(10000)

n_cells = 24
local_pool = 560
global_pool = 2000
n_synapses = 400

bin_size = 1

start_local = (global_pool - local_pool) / 2.0
# Draw local targets
local_targets = []
for x in range(n_cells):
    local_targets.append(np.random.choice(np.arange(start_local, start_local + local_pool), n_synapses, replace=False))

local_targets = np.array(local_targets)
plt.figure()
local_hist = plt.hist(local_targets.flatten(), bins = np.arange(-0.5,global_pool, bin_size))
plt.xlabel("# GC")
plt.ylabel("# Synapses")
plt.title("24 locally connecting cells")

# Draw global targets
global_targets = []
for x in range(n_cells):
    global_targets.append(np.random.choice(np.arange(0, global_pool), n_synapses, replace=False))

global_targets = np.array(global_targets)
plt.figure()
global_hist = plt.hist(global_targets.flatten(), bins = np.arange(-0.5,global_pool, bin_size))
plt.xlabel("# GC")
plt.ylabel("# Synapses")
plt.title("24 globally connecting cells")