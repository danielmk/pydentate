"""
Simulate and plot the intrinsic properties of a basket cell.

@author: danielmk
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
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
STEP_AMP_UP = 0.15  # in nA
STEP_AMP_DOWN = -0.2  # in nA

# Load compiled mechanisms
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Create Basket Cells
bc_one = basketcell.BasketCell()
bc_two = basketcell.BasketCell()

# Create Current Clamp for bc_one
cclamp_up = h.IClamp(bc_one.soma(0.5))
cclamp_up.delay = STEP_START
cclamp_up.dur = STEP_DUR
cclamp_up.amp = STEP_AMP_UP

# Create Current Clamp for bc_two
cclamp_down = h.IClamp(bc_two.soma(0.5))
cclamp_down.delay = STEP_START
cclamp_down.dur = STEP_DUR
cclamp_down.amp = STEP_AMP_DOWN

# Create Recording Vectors
bc_one_rec = h.Vector().record(bc_one.soma(0.5)._ref_v)
bc_two_rec = h.Vector().record(bc_two.soma(0.5)._ref_v)
clamp_up_rec = h.Vector().record(cclamp_up._ref_i)
clamp_down_rec = h.Vector().record(cclamp_down._ref_i)

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

ax[0].plot(time, bc_one_rec, color='k', label="Cell 1")
ax[0].plot(time, bc_two_rec, color='k', label="Cell 2")
ax[0].set_xlabel("Time (ms)")
ax[0].set_ylabel("Voltage (mV)")

ax[1].plot(time, clamp_up_rec, color='k', label="Cell 1")
ax[1].plot(time, clamp_down_rec, color='k', label="Cell 2")
ax[1].set_xlabel("Time (ms)")
ax[1].set_ylabel("Current (nA)")

plt.show()
