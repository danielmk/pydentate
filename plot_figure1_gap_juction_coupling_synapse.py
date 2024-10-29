"""
Simulate a spike in a granule cell that causes an EPSP in a basket cell that
has gap junction coupling with another basket cell.

Produces a plot showing the result.

@author: danielmk
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
import matplotlib.pyplot as plt
from pydentate import neuron_tools, basketcell, granulecell

# Constants
GAP_RESISTANCE = 600  # in Ohms
GAP_DELAY = 0  # in ms
DT = 0.1  # in ms
DURATION = 300  # in ms
WARMUP = 2000  # in ms
STEP_START = 100  # in ms
STEP_DUR = 30  # in ms
STEP_AMP = 0.2  # in nA

# Load compiled mechanisms
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Create Granule Cell and Basket Cells
gc = granulecell.GranuleCell()
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

# Create Synapse
syn = h.tmgsyn(bc_one.dendrites[1].secs[0](0.5))
syn.tau_1 = 8.7
syn.tau_facil = 500
syn.U = 0.1
syn.e = 0
syn.tau_rec = 0

# Create NetCon
THR = 10
DELAY = 0.8
WEIGHT = 0.025
netcon = h.NetCon(gc.soma(0.5)._ref_v, syn, THR, DELAY, WEIGHT, sec=gc.soma)

# Create Current Clamp
cclamp = h.IClamp(gc.soma(0.5))
cclamp.delay = STEP_START
cclamp.dur = STEP_DUR
cclamp.amp = STEP_AMP

# Create Recording Vectors
gc_rec = h.Vector().record(gc.soma(0.5)._ref_v)
bc_one_rec = h.Vector().record(bc_one.soma(0.5)._ref_v)
bc_two_rec = h.Vector().record(bc_two.soma(0.5)._ref_v)

# Run Neuron Simulator
neuron_tools.run_neuron_simulator(warmup=WARMUP, t_stop=DURATION, dt_sim=DT)

# Calculate Coupling Coefficient
cc_idx_one = np.array(bc_one_rec).argmax()
cc_idx_two = np.array(bc_two_rec).argmax()
bl_idx = int(50 / DT)
bl_one = np.mean(np.array(bc_one_rec)[:bl_idx])
bl_two = np.mean(np.array(bc_two_rec)[:bl_idx])
coupling_coefficient = (bc_two_rec[cc_idx_two] - bl_two) / (bc_one_rec[cc_idx_one] - bl_one)

# Plotting
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({'font.size': 22})
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']

time = np.linspace(0, DURATION, num=len(gc_rec))
fig, ax = plt.subplots(2, 1)

ax[0].plot(time, gc_rec, color=colors[2])
ax[0].set_xlabel("Time (ms)")
ax[0].set_ylabel("Voltage (mV)")

ax[1].plot(time, bc_one_rec, color=colors[0], label="Cell 1")
ax[1].plot(time, bc_two_rec, color=colors[1], label="Cell 2")
ax[1].set_xlabel("Time (ms)")
ax[1].set_ylabel("Voltage (mV)")
ax[1].legend()
ax[1].set_title(f"Coupling Coefficient: {coupling_coefficient:.4f}")

# Show the plot
plt.show()
