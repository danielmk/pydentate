# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 16:17:11 2018

@author: daniel
"""

from neuron import h, gui
import os

dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
            "C:\\Users\\daniel\\repos\\nrnmech.dll"]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

from granulecell import GranuleCell
from mossycell import MossyCell
from basketcell import BasketCell
from hippcell import HippCell
import matplotlib.pyplot as plt

gc = GranuleCell()
mc = MossyCell()
bc = BasketCell()
hc = HippCell()

gc_stim = gc._current_clamp_soma(0.8, 10, 5)
mc_stim = mc._current_clamp_soma(0.36, 500, 5)
bc_stim = bc._current_clamp_soma(0.5, 500, 50)
hc_stim = hc._current_clamp_soma(0.5, 500, 50)

gc_rec = gc._voltage_recording()
mc_rec = mc._voltage_recording()
bc_rec = bc._voltage_recording()
hc_rec = hc._voltage_recording()

time_vec = h.Vector()
time_vec.record(h._ref_t)
h.cvode.active(0)
dt = 0.1
h.steps_per_ms = 1.0/dt
h.tstop = 1000
h.finitialize(-70)
h.t = 0
h.secondorder = 0
h.dt = 10
while h.t < -100:
    h.fadvance()
h.secondorder = 2
h.t = 0
h.dt = 0.1

"""Setup run control for -100 to 1500"""
h.frecord_init() # Necessary after changing t to restart the vectors
while h.t < 600:
    h.fadvance()

plt.figure()
plt.plot(time_vec, gc_rec)

plt.figure()
plt.plot(time_vec, mc_rec)

plt.figure()
plt.plot(time_vec, bc_rec)

plt.figure()
plt.plot(time_vec, hc_rec)
