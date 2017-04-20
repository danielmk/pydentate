# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 11:16:23 2017

@author: DanielM
"""

from neuron import h, gui

class Network(object):

    def __init__(self, celltypes, cellnums):
        self.cells = {}
        for idx, x in enumerate(celltypes):
            self.cells[x] = []
            for y in range(cellnums[idx]):
                self.cells[x].append(x())

        self.pp_netstim = self.pp_stimulation

    def artificial_source(self, target_type, target_seg, Gmax, rise, decay, delay, stim):
        for cell in self.cells[target_type]:
            for dendrite in cell.dendrites:
                netstim = h.NetStim(0.5)
                netstim.interval = 100
                netstim.number = 1
                netstim.start = 5
                cell.netcon = h.NetCon(netstim.hocobjptr, dendrite[target_seg](0))
                #cell.pp_synapse.delay = delay

    def pp_stimulation(self, start, number = 1,  interval = 100):
        netstim = h.NetStim(0.5)
        netstim.start = start
        netstim.number = number
        netstim.interval = interval

        return netstim


