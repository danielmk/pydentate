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
import matplotlib.pyplot as plt

dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
             "C:\\Users\\daniel\\repos\\nrnmech.dll",
             "C:\\Users\\Daniel\\repos\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll"]
for x in dll_files:
    if os.path.isfile(x):
        dll_dir = x
print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll(dll_dir)

current_steps = np.arange(0, 0.6, 0.025)
gc_voltages = []
mc_voltages = []
bc_voltages = []
hc_voltages = []

gc_apcs = []
mc_apcs = []
bc_apcs = []
hc_apcs = []

dt = 0.1
for x in current_steps:

    gc = GranuleCell()
    mc = MossyCell()
    bc = BasketCell()
    hc = HippCell()

    for sec in gc.all_secs:
        sec.ek = sec.ek+0
    for sec in mc.all_secs:
        sec.ek = sec.ek+0
    for sec in bc.all_secs:
        sec.ek = sec.ek+0
    for sec in hc.all_secs:
        sec.ek = sec.ek+0

    gc_stim = gc._current_clamp_soma(x, 500, 50)
    mc_stim = mc._current_clamp_soma(x, 500, 50)
    bc_stim = bc._current_clamp_soma(x, 500, 50)
    hc_stim = hc._current_clamp_soma(x, 500, 50)

    gc_rec = gc._voltage_recording()
    mc_rec = mc._voltage_recording()
    bc_rec = bc._voltage_recording()
    hc_rec = hc._voltage_recording()

    gc_apc_curr = h.APCount(gc.soma(0.5))
    mc_apc_curr = h.APCount(mc.soma(0.5))
    bc_apc_curr = h.APCount(bc.soma(0.5))
    hc_apc_curr = h.APCount(hc.soma(0.5))

    time_vec = h.Vector()
    time_vec.record(h._ref_t)

    h.cvode.active(0)
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

    gc_apcs.append(gc_apc_curr.n / 0.5)
    mc_apcs.append(mc_apc_curr.n / 0.5)
    bc_apcs.append(bc_apc_curr.n / 0.5)
    hc_apcs.append(hc_apc_curr.n / 0.5)

plt.figure()
plt.plot(current_steps*1000, gc_apcs, marker = 'o')
plt.title("GC F-I Curve")
plt.xlabel("Current Injection (pA)")
plt.ylabel("Frequency (Hz)")
plt.figure()
plt.plot(current_steps*1000, mc_apcs, marker = 'o')
plt.title("MC F-I Curve")
plt.xlabel("Current Injection (pA)")
plt.ylabel("Frequency (Hz)")
plt.figure()
plt.plot(current_steps*1000, bc_apcs, marker = 'o')
plt.title("BC F-I Curve")
plt.xlabel("Current Injection (pA)")
plt.ylabel("Frequency (Hz)")
plt.figure()
plt.plot(current_steps*1000, hc_apcs, marker = 'o')
plt.title("HC F-I Curve")
plt.xlabel("Current Injection (pA)")
plt.ylabel("Frequency (Hz)")
