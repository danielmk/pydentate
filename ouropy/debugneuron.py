# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 11:16:20 2017

@author: DanielM
"""

from ouropy.genneuron import GenNeuron
import ouropy.parameters as params
from neuron import h, gui
import matplotlib.pyplot as plt

class DebugNeuron(GenNeuron):
    """Creates a neuron for debugging purposes"""
    
    def __init__(self):
        self.mk_soma('soma', 16.8,16.8)
        self.mk_dendrite(4, dend_name = 'dend_1', sec_names = ['gcld', 'proxd', 'midd', 'dd'], diam = [3,3,3,3], L = [50,150,150,150])
        self.mk_dendrite(4, dend_name = 'dend_2', sec_names = ['gcld', 'proxd', 'midd', 'dd'], diam = [3,3,3,3], L = [50,150,150,150])

        

if __name__ == '__main__':
    mydn = DebugNeuron()
    GCparams = params.read_parameters("C:\\Users\\DanielM\\Repos\\pyDentate\\granulecellparams.txt")
    mydn.insert_mechs(GCparams)
    mydn._current_clamp_soma(delay = 200)
    voltage, time = mydn._voltage_recording()
    h.tstop = 800.0
    h.run()
    plt.plot(time, voltage)