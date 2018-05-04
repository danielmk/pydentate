# -*- coding: utf-8 -*-
"""
Created on Wed May 02 15:21:55 2018

@author: DanielM
"""

from neuron import h
import ouropy
import numpy as np
import net_tuned
import net_nonfacilitating
import net_global
#import net_tuned_10ECInputs
import os
import matplotlib.pyplot as plt

n_cells = 2000
n_pps = 400

np.random.seed(10000)
# Original from Yim
#np.random.seed(10000)
#temporal_patterns = np.random.poisson(10, (400, 3)).cumsum(axis = 1)

# Generate the PP -> GC mapping so that each GC receives inputs from 20/400
# randomly chosen PP inputs
innervation_pattern_gc = np.array([np.random.choice(n_pps,20, replace = False) for x in range(n_cells)])
innervation_pattern_gc = innervation_pattern_gc.swapaxes(0,1)

PP_to_GCs = []
for x in range(0,n_pps):
    PP_to_GCs.append(np.argwhere(innervation_pattern_gc == x)[:,1])

PP_to_GCs = np.array(PP_to_GCs)

average_convergence_list = []
for x in range(n_pps-24):
    hist = np.histogram(np.hstack(PP_to_GCs[0+x:24+x].flat), bins=np.arange(-0.5,n_cells+0.5,1))
    average_convergence_list.append(hist[0].mean())

average_convergence = np.array(average_convergence_list).mean()
std_convergence = np.array(average_convergence_list).std()

"""hist = plt.hist(np.hstack(PP_to_GCs[0:6].flat), bins=np.arange(-0.5,n_cells+0.5,1))
average_convergence_list.append(hist[0].mean)"""