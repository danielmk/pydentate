# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 16:17:11 2018

@author: daniel
"""

from neuron import h, gui
import os
import numpy as np
from granulecell import GranuleCell
from mossycell import MossyCell
from basketcell import BasketCell
from hippcell import HippCell

dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
             "C:\\Users\\daniel\\repos\\nrnmech.dll",
             "C:\\Users\\Daniel\\repos\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll"]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

current_steps = np.arange(-0.3, 1, 0.05)
gc_voltages = []
mc_voltages = []
bc_voltages = []
hc_voltages = []

for x in current_steps:

    gc = GranuleCell()
    mc = MossyCell()
    bc = BasketCell()
    hc = HippCell()

    gc_stim = gc._current_clamp_soma(x, 500, 50)
    mc_stim = mc._current_clamp_soma(x, 500, 50)
    bc_stim = bc._current_clamp_soma(x, 500, 50)
    hc_stim = hc._current_clamp_soma(x, 500, 50)

    gc_rec = gc._voltage_recording()
    mc_rec = mc._voltage_recording()
    bc_rec = bc._voltage_recording()
    hc_rec = hc._voltage_recording()

    time_vec = h.Vector()
    time_vec.record(h._ref_t)

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
    h.frecord_init()  # Necessary after changing t to restart the vectors

    while h.t < 700:
        h.fadvance()

    gc_voltages.append(np.array(gc_rec))
    mc_voltages.append(np.array(mc_rec))
    bc_voltages.append(np.array(bc_rec))
    hc_voltages.append(np.array(hc_rec))
