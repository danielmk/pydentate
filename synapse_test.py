# -*- coding: utf-8 -*-
"""
Created on Mon Mar 05 13:41:23 2018

@author: DanielM
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace 
import numpy as np
import ouropy
from input_generator import inhom_poiss
import os
import argparse
import time
from analysis_main import time_stamps_to_signal
import pdb
from ouropy.gennetwork import GenNetwork, ConnDivergent, ImplicitConvergentTmgsynConnectionExpProb
from granulecell import GranuleCell
import matplotlib.pyplot as plt

# Where to search for nrnmech.dll file. Must be adjusted for your machine.
dll_files = [("/home/daniel/repos/pyDentate/mechs_7-6_linux/x86_64/.libs/libnrnmech.so"),
             ("C:\\Users\\Daniel\\repos\\pyDentate\\mechs_7-6_win\\nrnmech.dll")]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
h.nrn_load_dll(dll_dir)

seed=1000

nw_ppdynamicstuned = GenNetwork()
nw_ppdynamicstuned.mk_population(GranuleCell, 1)

nw_ppdynamicsoff = GenNetwork()
nw_ppdynamicsoff.mk_population(GranuleCell, 1)

ctuned = 'g'
coff = 'b'
t_patterns = np.arange(50,1001,200)
t_patterns = t_patterns[np.newaxis,:]
W_pp_gc = 1e-3

ImplicitConvergentTmgsynConnectionExpProb(nw_ppdynamicstuned.populations[0], t_patterns, 'midd', 1,
             1.5, 10, 2.22590184e+13, 0.1, 3.12395254e+02, 8.15226996e-01, 0, W_pp_gc, rec_cond=True, n_inputs=1)

ImplicitConvergentTmgsynConnectionExpProb(nw_ppdynamicsoff.populations[0], t_patterns, 'midd', 1,
             1.5, 10, 0, 0.1, 0, 0, 0, W_pp_gc, rec_cond=True, n_inputs=1)

"""Initialization for -2000 to -100"""
h.cvode.active(0)
dt = 0.1
h.steps_per_ms = 1.0/dt
h.finitialize(-60)
h.t = -2000
h.secondorder = 0
h.dt = 10
while h.t < -100:
    h.fadvance()

h.secondorder = 2
h.t = 0
h.dt = 0.1

"""Setup run control for -100 to 1500"""
h.frecord_init()  # Necessary after changing t to restart the vectors

while h.t < 1100:
    h.fadvance()

"""Plotting"""
x = np.arange(0,1100.2,0.1)
plt.plot(x, np.array(nw_ppdynamicstuned.populations[0].connections[0].conductances[0][1]), color=ctuned)
plt.plot(x, np.array(nw_ppdynamicsoff.populations[0].connections[0].conductances[0][1]), color=coff)
plt.legend(("ppdynamicstuned", "ppdynamicsoff"))
plt.ylabel("synaptic conductance (uS)")
plt.xlabel("time (ms)")