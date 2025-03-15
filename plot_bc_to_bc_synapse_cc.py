"""
Simulate a BC-BC synapse and measure postsynaptic voltage deflection. NOT DONE!

@author: danielmk
"""

from neuron import h, gui  # gui necessary for some parameters to h namespace
import numpy as np
from pydentate import neuron_tools
from pydentate import  basketcell
import matplotlib.pyplot as plt

# Parameters
dt = 0.01  # In ms
duration = 300  # In ms
warmup = 2000  # In ms

step_start = 100  # In ms
step_dur = 20  # In ms
step_amp = 0.2 # In nA

neuron_tools.load_compiled_mechanisms(path=r'C:\Users\Daniel\repos\pydentate\mechs\nrnmech.dll')

bc_one = basketcell.BasketCell()
bc_two = basketcell.BasketCell()

# Create Synapse
syn = h.tmgsyn(bc_two.dendrites[1].secs[0](0.5))
syn.tau_1 = 1.8
syn.tau_facil = 0
syn.U = 1
syn.e = -70
syn.tau_rec = 0
thr = 10
delay = 0.8
weight = 7.6e-3
# weight = 15.2e-3
netcon = h.NetCon(bc_one.soma(0.5)._ref_v, syn, thr, delay, weight, sec=bc_one.soma)

cclamp = h.IClamp(bc_one.soma(0.5))

cclamp.delay = step_start
cclamp.dur = step_dur
cclamp.amp = step_amp


cclamp_record = h.Vector()
cclamp_record.record(bc_two.soma(0.5)._ref_v)

bc_one_rec = h.Vector()
bc_one_rec.record(bc_one.soma(0.5)._ref_v)

neuron_tools.run_neuron_simulator(warmup=warmup, t_stop=duration, dt_sim=dt)


# vclamp_record


# Plotting
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams.update({'font.size': 22})
colors = ['#C70019', '#0D6B9A', '#EE9A20', '#6389A5', '#EA521C', '#8A963F']

fig, ax = plt.subplots(1, 1)

ax.plot(t, bc_one_rec, color=colors[0])
ax.plot(t, cclamp_record, color=colors[1])
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Voltage (mV)")
ax.set_xlim((100, 150))


