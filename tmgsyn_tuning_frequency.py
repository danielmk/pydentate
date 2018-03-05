# -*- coding: utf-8 -*-
"""
Created on Thu Jan 04 15:26:50 2018

@author: DanielM
"""

from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt
from pyDentate.mossycell_cat import MossyCell
from pyDentate.basketcell import BasketCell

h.nrn_load_dll("C:\\Users\\DanielM\\Repos\\models_dentate\\dentate_gyrus_Santhakumar2005_and_Yim_patterns\\dentategyrusnet2005\\nrnmech.dll")

stim_periods = [1000, 100, 33, 20]

for period in stim_periods:
    """Setup stimulation pattern"""
    t_pattern = np.arange(100, 100+10*period, period)
    vecstim = h.VecStim()
    pattern_vec = h.Vector(t_pattern)
    vecstim.play(pattern_vec)

    """Setup tmgsyn"""
    mc_tmgsyn = MossyCell()

    mc_tmgsyn_syn = h.tmgsyn(mc_tmgsyn.all_secs[1](0.5))
    mc_tmgsyn_syn.e = 0
    mc_tmgsyn_syn.tau_facil = 500 # This parameter gives the frequency dependence of facilitation
    mc_tmgsyn_syn.tau_1 = 6.2
    mc_tmgsyn_syn.tau_rec = 0 # ???
    mc_tmgsyn_syn.U = 0.1
    #mc_tmgsyn_syn.u0 = 0.04
    mc_tmgsyn_netcon = h.NetCon(vecstim, mc_tmgsyn_syn)
    mc_tmgsyn_netcon.weight[0] = 0.2*10**(-3) * 10

    mc_tmgsyn_rec_g = h.Vector()
    mc_tmgsyn_rec_g.record(mc_tmgsyn_syn._ref_g)

    mc_tmgsyn._SEClamp(dur1=t_pattern[9]+500, amp1=-70, rs=0.001)
    rec_curr = h.Vector()
    rec_curr.record(mc_tmgsyn.vclamp._ref_i)

    #mc_tmgsyn_syn.U -> facilitation coefficient
    #mc_tmgsyn_syn.tau_facil -> time decay of facilitation
    #mc_tmgsyn_syn.tau_1 -> decay constant of synaptic conductance
    #mc_tmgsyn_syn.tau_rec -> ???

    dt = 0.01
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
    h.dt = 0.01

    """Setup run control for -100 to 1500"""
    h.frecord_init() # Necessary after changing t to restart the vectors
    while h.t < t_pattern[9]+500:
        h.fadvance()

    plt.figure()
    plt.plot(rec_curr)
    #plt.plot(mc_exp2syn_rec_v, color = 'g')    



