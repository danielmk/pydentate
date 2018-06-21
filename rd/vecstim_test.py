# -*- coding: utf-8 -*-
"""
Created on Mon Dec 04 15:32:13 2017

@author: DanielM
"""

from neuron import h, gui
import matplotlib.pyplot as plt

h.nrn_load_dll("C:\Users\DanielM\Repos\models_dentate\dentate_gyrus_Santhakumar2005_and_Yim_patterns\dentategyrusnet2005\\nrnmech.dll")

soma = h.Section(name='soma')
ampa = h.Exp2Syn(0.5, sec=soma)
cond_vec = h.Vector()
cond_vec.record(ampa._ref_g)

vecStim = h.VecStim()
vec = h.Vector([1, 30])
vecStim.play(vec)

netCon = h.NetCon(vecStim, ampa)
netCon.weight[0] = 1
time_vec = h.Vector()
time_vec.record(h._ref_t)

h.tstop = 100
h.run()