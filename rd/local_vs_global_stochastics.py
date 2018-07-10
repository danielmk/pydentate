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

# Draw local targets
local_targets = []
for x in range(n_cells):
    local_targets.append(np.random.choice(np.arange(local_pool), 400, replace=False))
    
local_targets = np.array(local_targets)
