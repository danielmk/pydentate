# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 16:17:11 2018

@author: daniel
"""


import matplotlib.pyplot as plt
import numpy as np
from neuron import h

h.load_file("stdrun.hoc")

from pydentate import BasketCell, GranuleCell, HippCell, MossyCell

# dll_files = ["C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll",
#              "C:\\Users\\daniel\\repos\\nrnmech.dll",
#              "C:\\Users\\Daniel\\repos\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll"]
# for x in dll_files:
#     if os.path.isfile(x):
#         dll_dir = x
# print("DLL loaded from: " + str(dll_dir))
h.nrn_load_dll("./pydentate/x86_64/.libs/libnrnmech.so")

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
    h.steps_per_ms = 1.0 / dt
    h.tstop = 1500
    h.finitialize(-60)
    h.t = -2000
    h.secondorder = 0
    h.dt = 10
    while h.t < -100:
        h.fadvance()
        # print(h.t)

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

fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(8.27, 11.69))

axes[0].plot(current_steps * 1000, gc_apcs, marker="o")
axes[0].set_title("GC F-I Curve")
axes[0].set_xlabel("Current Injection (pA)")
axes[0].set_ylabel("Frequency (Hz)")

axes[1].plot(current_steps * 1000, mc_apcs, marker="o")
axes[1].set_title("MC F-I Curve")
axes[1].set_xlabel("Current Injection (pA)")
axes[1].set_ylabel("Frequency (Hz)")

axes[2].plot(current_steps * 1000, bc_apcs, marker="o")
axes[2].set_title("BC F-I Curve")
axes[2].set_xlabel("Current Injection (pA)")
axes[2].set_ylabel("Frequency (Hz)")

axes[3].plot(current_steps * 1000, hc_apcs, marker="o")
axes[3].set_title("HC F-I Curve")
axes[3].set_xlabel("Current Injection (pA)")
axes[3].set_ylabel("Frequency (Hz)")

fig.tight_layout()
fig.savefig("intrinsic_properties.png")
