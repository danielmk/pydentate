# -*- coding: utf-8 -*-
"""
Created on Thu Jan 04 15:26:50 2018

@author: DanielM
"""

from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from pyDentate.mossycell_cat import MossyCell

h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")

"""Setup stimulation pattern"""
t_pattern = np.arange(50,350,33)
vecstim = h.VecStim()
pattern_vec = h.Vector(t_pattern)
vecstim.play(pattern_vec)

"""Setup tmgsyn"""
mc_tmgsyn = MossyCell()

mc_tmgsyn_syn = h.tmgsyn(mc_tmgsyn.all_secs[1](0.5))
#curr_syn.tau1 = tau1
#curr_syn.e = 0
#curr_syn.tau_facil = 33 # This parameter gives the frequency dependence of facilitation
#curr_syn.tau_1 = 1
mc_tmgsyn_netcon = h.NetCon(vecstim, mc_tmgsyn_syn)
mc_tmgsyn_netcon.weight[0] = 0.05

mc_tmgsyn_rec_g = h.Vector()
mc_tmgsyn_rec_g.record(mc_tmgsyn_syn._ref_g)
mc_tmgsyn_rec_v = h.Vector()
mc_tmgsyn_rec_v.record(mc_tmgsyn.soma(0.5)._ref_v)

"""Setup Exp2Syn"""
mc_exp2syn = MossyCell()

mc_exp2syn_syn = h.Exp2Syn(mc_exp2syn.all_secs[1](0.5))
mc_exp2syn_netcon = h.NetCon(vecstim, mc_exp2syn_syn)
mc_exp2syn_netcon.weight[0] = 0.05

mc_exp2syn_rec_g = h.Vector()
mc_exp2syn_rec_g.record(mc_exp2syn_syn._ref_g)
mc_exp2syn_rec_v = h.Vector()
mc_exp2syn_rec_v.record(mc_exp2syn.soma(0.5)._ref_v)


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