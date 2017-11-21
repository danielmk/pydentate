# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 09:49:04 2017

@author: DanielM
"""

from granulecell import GranuleCell
from neuron import h, gui
import matplotlib.pyplot as plt
import numpy as np
import scalebars as sb

h.nrn_load_dll("C:\Users\DanielM\Repos\models_dentate\dentate_gyrus_Santhakumar2005_corrected\\nrnmech.dll")

"""Setup cell and measurements"""
gc = GranuleCell()

soma_v = h.Vector()
soma_cai = h.Vector()
soma_ik_cagk = h.Vector()
soma_ik_cagk = h.Vector()

t = h.Vector()
soma_v.record(gc.soma(0.5)._ref_v)
soma_cai.record(gc.soma(0.5)._ref_ncai)
soma_ik_cagk.record(gc.soma(0.5)._ref_ik_cagk)

stim = h.IClamp(gc.soma(0.5))
stim.delay = 300
stim.amp = 0.2
stim.dur = 500
t.record(h._ref_t)

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
while h.t < 1500:
    h.fadvance()

"""Setup plots"""
plt.figure()
plt.plot(t, soma_v)
plt.figure()
plt.plot(t,soma_cai)
plt.figure()
plt.plot(t,soma_ik_cagk)


v_array = soma_v.as_numpy()
