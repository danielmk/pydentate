# -*- coding: utf-8 -*-
"""
Created on Thu Jan 04 15:26:50 2018

@author: DanielM
"""

from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from pyDentate.granulecell import GranuleCell
from pyDentate.mossycell_cat import MossyCell

h.nrn_load_dll("C:\\Users\\daniel\\repos\\nrnmech.dll")

post = h.Section()

t_pattern = np.arange(50,350,33)

mc = MossyCell()

vecstim = h.VecStim()
pattern_vec = h.Vector(t_pattern)
vecstim.play(pattern_vec)

curr_syn = h.tmgsyn(mc.all_secs[1](0.5))
#curr_syn.tau1 = tau1
#curr_syn.tau2 = tau2
curr_syn.e = 0
curr_syn.tau_facil = 30 # This parameter gives the frequency dependence of facilitation
curr_netcon = h.NetCon(vecstim, curr_syn)
curr_netcon.weight[0] = 0.05

rec_g = h.Vector()
rec_g.record(curr_syn._ref_g)
rec_v = h.Vector()
rec_v.record(mc.soma(0.5)._ref_v)

dt = 0.1
h.steps_per_ms = 1.0/dt
h.tstop = 1500
h.finitialize(-60)
h.t = -2000
h.secondorder = 0
h.dt = 10
while h.t < -100:
    h.fadvance()
    #print(h.t)
h.secondorder = 2
h.t = 0
h.dt = 0.1

"""Setup run control for -100 to 1500"""
h.frecord_init() # Necessary after changing t to restart the vectors
while h.t < 500:
    h.fadvance()