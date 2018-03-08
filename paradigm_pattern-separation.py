# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h
import ouropy
import numpy as np
import net_tuned
import time
# Office PC
#h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")
#Home PC
h.nrn_load_dll("C:\\Users\\daniel\\repos\\nrnmech.dll")
np.random.seed(10000)
# Generate temporal patterns for the 100 PP inputs
temporal_patterns = np.random.poisson(10, (400, 3)).cumsum(axis = 1)

# Generate the PP -> GC mapping so that each GC receives inputs from 20/400
# randomly chosen PP inputs
innervation_pattern_gc = np.array([np.random.choice(400,20, replace = False) for x in range(2000)])
innervation_pattern_gc = innervation_pattern_gc.swapaxes(0,1)

PP_to_GCs = []
for x in range(0,400):
    PP_to_GCs.append(np.argwhere(innervation_pattern_gc == x)[:,1])

PP_to_GCs = np.array(PP_to_GCs)

innervation_pattern_bc = np.array([np.random.choice(400,20, replace = False) for x in range(24)])
innervation_pattern_bc = innervation_pattern_bc.swapaxes(0,1)

PP_to_BCs = []
for x in range(0,400):
    PP_to_BCs.append(np.argwhere(innervation_pattern_bc == x)[:,1])

PP_to_BCs = np.array(PP_to_BCs)
all_targets = np.array([y for x in PP_to_GCs for y in x])

nw = net_tuned.TunedNetwork(10000, temporal_patterns[0:6], PP_to_GCs[0:6], PP_to_BCs[0:6], sprouting=0)

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

while h.t < 200:
    h.fadvance()