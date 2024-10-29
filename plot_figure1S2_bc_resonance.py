"""
Simulate and plot the intrinsic properties of a basket cell.

@author: danielmk
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
import matplotlib.pyplot as plt
from pydentate import neuron_tools, basketcell
from scipy.signal import chirp

# Constants
GAP_RESISTANCE = 600  # in Ohms
GAP_DELAY = 0  # in ms
DT = 0.1  # in ms
DURATION = 5000  # in ms
WARMUP = 2000  # in ms
STEP_START = 100  # in ms
STEP_DUR = 500  # in ms
STEP_AMP_UP = 0.15  # in nA
STEP_AMP_DOWN = -0.2  # in nA

# Load compiled mechanisms
neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

# Create Basket Cells
bc = basketcell.BasketCell()

# Create Current Clamp Vector of Increasing Frequency
input_scale = 0.1
time = np.arange(0, DURATION + 1 * DT, DT)
chirp_input = chirp(time/1000, 1, DURATION/1000, 100) * input_scale
chirp_vector = h.Vector(chirp_input)
time_vector = h.Vector(time)

# Create Current Clamp for bc_one
# cclamp = h.IClamp(bc.soma(0.5))
cclamp = h.IClamp(bc.soma(0.5))
cclamp.delay = 0
cclamp.dur = DURATION

chirp_vector.play(cclamp._ref_amp, time_vector)

# Create Recording Vectors
bc_rec = h.Vector().record(bc.soma(0.5)._ref_v)
clamp_rec = h.Vector().record(cclamp._ref_i)

# Run Neuron Simulator
neuron_tools.run_neuron_simulator(warmup=WARMUP, t_stop=DURATION, dt_sim=DT)

# Plotting
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({'font.size': 22})
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']

fig, ax = plt.subplots(2, 1)

ax[0].plot(time, bc_rec, color='k', label="Cell 1")
ax[0].set_xlabel("Time (ms)")
ax[0].set_ylabel("Voltage (mV)")

ax[1].plot(time, clamp_rec, color='k', label="Cell 1")
ax[1].set_xlabel("Time (ms)")
ax[1].set_ylabel("Current (nA)")

plt.show()
