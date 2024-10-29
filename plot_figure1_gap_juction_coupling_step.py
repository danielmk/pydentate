"""
Simulate a step current injection into a basket cell that has gap junction
coupling with another basket cell.

Produces a plot showing the result.

@author: danielmk
"""

from neuron import h, gui  # gui necessary for some parameters in h namespace
import numpy as np
import matplotlib.pyplot as plt
from pydentate import neuron_tools, basketcell

# Constants
GAP_RESISTANCE = 600  # in Ohms
GAP_DELAY = 0  # in ms
DT = 0.1  # in ms
DURATION = 800  # in ms
WARMUP = 2000  # in ms
STEP_START = 100  # in ms
STEP_DUR = 500  # in ms
STEP_AMP = 0.1  # in nA

# Load compiled mechanisms
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Create Basket Cells
bc_one = basketcell.BasketCell()
bc_two = basketcell.BasketCell()

# Create Gap Junctions
bc_one_gap = h.gap(bc_one.dendrites[1].secs[0](0.5))
bc_two_gap = h.gap(bc_two.dendrites[1].secs[0](0.5))

# Set Gap Junction Parameters
bc_one_gap.r = GAP_RESISTANCE
bc_two_gap.r = GAP_RESISTANCE
bc_one_gap.delay = GAP_DELAY
bc_two_gap.delay = GAP_DELAY

# Set up Gap Junction Pointers
h.setpointer(bc_one.dendrites[1].secs[0](0.5)._ref_v, "v_pair", bc_two_gap)
h.setpointer(bc_two.dendrites[1].secs[0](0.5)._ref_v, "v_pair", bc_one_gap)

# Create Current Clamp
cclamp = h.IClamp(bc_one.soma(0.5))
cclamp.delay = STEP_START
cclamp.dur = STEP_DUR
cclamp.amp = STEP_AMP

# Create Recording Vectors
bc_one_rec = h.Vector().record(bc_one.soma(0.5)._ref_v)
bc_two_rec = h.Vector().record(bc_two.soma(0.5)._ref_v)
clamp_rec = h.Vector().record(cclamp._ref_i)

# Run Neuron Simulator
neuron_tools.run_neuron_simulator(warmup=WARMUP, t_stop=DURATION, dt_sim=DT)

# Calculate Coupling Coefficient
cc_idx = int(500 / DT)
bl_idx = int(50 / DT)
bl_one = np.mean(np.array(bc_one_rec)[0:bl_idx])
bl_two = np.mean(np.array(bc_two_rec)[0:bl_idx])
coupling_coefficient = (bc_two_rec[cc_idx] - bl_two) / (bc_one_rec[cc_idx] - bl_one)

# Plotting
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({'font.size': 22})
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']

fig, ax = plt.subplots(2, 1)
time = np.arange(0, DURATION + 2 * DT, DT)

ax[0].plot(time, bc_one_rec, color=colors[0], label="Cell 1")
ax[0].plot(time, bc_two_rec, color=colors[1], label="Cell 2")
ax[0].set_xlabel("Time (ms)")
ax[0].set_ylabel("Voltage (mV)")
ax[0].set_title(f"Coupling Coefficient: {coupling_coefficient:.4f}")

ax[1].plot(time, clamp_rec, color=colors[0])
ax[1].set_xlabel("Time (ms)")
ax[1].set_ylabel("Current (nA)")

# Show the plot
plt.show()

