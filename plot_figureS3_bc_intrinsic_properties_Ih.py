"""
Simulate and plot the intrinsic properties of a basket cell.

@author: danielmk
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
from pydentate import neuron_tools
from pydentate import  basketcell
import matplotlib.pyplot as plt

# Define Constants
gap_resistance = 6e2
gap_delay = 0
dt = 0.1  # In ms
duration = 2000  # In ms
warmup = 2000  # In ms

step_start = 100  # In ms
step_dur = 1000  # In ms
step_amp = 0.15  # In nA

# Load compiled mechanisms
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Create Basket Cells
bc_one = basketcell.BasketCell_Ih()
bc_two = basketcell.BasketCell_Ih()

"""Insert Ih and adjust intrinsic properties"""
bc_one.soma.insert('Ih')
bc_two.soma.insert('Ih')
bc_one.soma(0.5).Ih.gkhbar = 0.01762
bc_two.soma(0.5).Ih.gkhbar = 0.01762


# Create Current Clamp
cclamp_up = h.IClamp(bc_one.soma(0.5))
cclamp_up.delay = step_start
cclamp_up.dur = step_dur
cclamp_up.amp = step_amp

cclamp_down = h.IClamp(bc_two.soma(0.5))
cclamp_down.delay = step_start
cclamp_down.dur = step_dur
cclamp_down.amp = -0.1


# Create Recording Vectors
bc_one_rec = h.Vector()
bc_one_rec.record(bc_one.soma(0.5)._ref_v)
bc_two_rec = h.Vector()
bc_two_rec.record(bc_two.soma(0.5)._ref_v)

clamp_up_rec = h.Vector()
clamp_up_rec.record(cclamp_up._ref_i)
clamp_down_rec = h.Vector()
clamp_down_rec.record(cclamp_down._ref_i)

# Run Neuron Simulator
neuron_tools.run_neuron_simulator(warmup=warmup, t_stop=duration, dt_sim=dt)

"""CALCULATE THE SAG RATIO"""
bl = np.array(bc_two_rec)[0:100].mean(0)
ss = np.array(bc_two_rec)[10000:10050].mean(0)
sag = np.array(bc_two_rec).min()

sag_ratio = (bl - ss) / (bl - sag)

print(f"Sag Ratio: {sag_ratio}")


# Plotting
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({'font.size': 22})
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']

fig, ax = plt.subplots(3, 1)
t = np.arange(0, 2000.2, 0.1)
ax[0].plot(t, bc_one_rec, color='k', label="Cell 1")
ax[0].set_xlabel("Time (ms)")
ax[0].set_ylabel("Voltage (mV)")

ax[1].plot(t, bc_two_rec, color='k', label="Cell 2")
ax[1].set_xlabel("Time (ms)")
ax[1].set_ylabel("Voltage (mV)")

ax[2].plot(t, clamp_up_rec, color='k', label="Cell 1")
ax[2].plot(t, clamp_down_rec, color='k', label="Cell 2")
ax[2].set_xlabel("Time (ms)")
ax[2].set_ylabel("Current (nA)")
