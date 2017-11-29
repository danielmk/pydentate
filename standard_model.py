# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 13:01:38 2017

@author: DanielM
"""

from neuron import h, gui
from standardnetwork import StandardNetwork
import matplotlib.pyplot as plt
import numpy as np

nw = StandardNetwork(seed = 10000, sprouting = 0)

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
    #print(h.t)
h.secondorder = 2
h.t = 0
h.dt = 0.1

"""Setup run control for -100 to 1500"""
h.frecord_init() # Necessary after changing t to restart the vectors
while h.t < 300:
    h.fadvance()