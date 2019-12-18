# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 09:30:06 2019

@author: Daniel
"""

from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt

# locate nrnmech.dll and load it into h
dll_dir = ("C:\\Users\\Daniel\\Dropbox\\MEC Project\\BarisProject\\tmgexp2syn"
           "\\nrnmech.dll")
h.nrn_load_dll(dll_dir)

# Set some parameters
tau_1 = 2
tau_2 = 20
gmax = 1
pas_e = -65
syn_e = 0
tau_facil = 500
tau_rec = 0
stim_interval = 200
stim_number = 10
stim_noise = 0
U = 0.1
weight=1

# Create the sections where the synapses will be attached
exp2syn_sec = h.Section()
tmgsyn_sec = h.Section()
tmgexp2syn_sec = h.Section()
secs = [exp2syn_sec, tmgsyn_sec, tmgexp2syn_sec]
for sec in secs:
    sec.insert('pas')
    sec(0.5).pas.e = pas_e

# Create the synapses
exp2syn_syn = h.Exp2Syn(exp2syn_sec(0.5))
exp2syn_syn.tau1 = tau_1
exp2syn_syn.tau2 = tau_2
exp2syn_syn.e = syn_e
exp2syn_syn.g = gmax

tmgsyn_syn = h.tmgsyn(tmgsyn_sec(0.5))
tmgsyn_syn.tau_1 = tau_2  # Note that tau_2 is the decay tau at script level!
tmgsyn_syn.tau_facil = tau_facil
tmgsyn_syn.tau_rec = tau_rec
tmgsyn_syn.e = syn_e
#tmgsyn_syn.g = gmax
tmgsyn_syn.U = U

tmgexp2syn_syn = h.tmgexp2syn(tmgexp2syn_sec(0.5))
tmgexp2syn_syn.tau_1 = tau_1
tmgexp2syn_syn.tau_2 = tau_2
tmgexp2syn_syn.tau_facil = tau_facil
tmgexp2syn_syn.tau_rec = tau_rec
tmgexp2syn_syn.e = syn_e
#tmgexp2syn_syn.g = gmax
tmgexp2syn_syn.U = U

syns = [exp2syn_syn, tmgsyn_syn, tmgexp2syn_syn]

# Create the netstim object that feeds events into the synapses
netstim = h.NetStim()
netstim.start = 50
netstim.interval = stim_interval
netstim.number = stim_number
netstim.noise = stim_noise

# Connect the netstim with the synapses
exp2syn_netcon = h.NetCon(netstim, exp2syn_syn)
exp2syn_netcon.weight[0] = weight
tmgsyn_netcon = h.NetCon(netstim, tmgsyn_syn)
tmgsyn_netcon.weight[0] = weight
tmgexp2syn_netcon = h.NetCon(netstim, tmgexp2syn_syn)
tmgexp2syn_netcon.weight[0] = weight

# Setup the conductance recordings
exp2syn_rec = h.Vector()
exp2syn_rec.record(exp2syn_syn._ref_g)

tmgsyn_rec = h.Vector()
tmgsyn_rec.record(tmgsyn_syn._ref_g)

tmgexp2syn_rec = h.Vector()
tmgexp2syn_rec.record(tmgexp2syn_syn._ref_g)

# Run the simulation
h.cvode.active(0)
dt = 0.1
h.steps_per_ms = 1.0/dt
h.finitialize(-65)
h.t = -2000
h.secondorder = 0
h.dt = 10
while h.t < -100:
    h.fadvance()

h.secondorder = 2
h.t = 0
h.dt = 0.1

"""Setup run control for -100 to 1500"""
h.frecord_init()  # Necessary after changing t to restart the vectors

while h.t < 500:
    h.fadvance()

# Plotting
exp2syn_rec = np.array(exp2syn_rec)
tmgsyn_rec = np.array(tmgsyn_rec)
tmgexp2syn_rec = np.array(tmgexp2syn_rec)

plt.figure()
plt.plot(exp2syn_rec)
plt.title("exp2syn")
plt.ylim((0,0.6))

plt.figure()
plt.plot(tmgsyn_rec)
plt.ylim((0,0.6))
plt.plot(tmgexp2syn_rec)
plt.legend(("g_tmgsyn", "g_tmgexp2syn"))
plt.ylim((0,0.6))

plt.figure()
plt.plot(exp2syn_rec/exp2syn_rec.max())