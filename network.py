# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 11:16:23 2017

@author: DanielM
"""

from neuron import h, gui
from granulecell import GranuleCell

class Network(object):

    def __init__(self, celltypes, cellnums):
        self.cells = {}
        for idx, x in enumerate(celltypes):
            self.cells[x] = []
            for y in range(cellnums[idx]):
                self.cells[x].append(x())

        #self.pp_netstim = self.pp_stimulation

    def pp_stimulation(self, start, number = 1,  interval = 100):
        netstim = h.NetStim(0.5)
        netstim.start = start
        netstim.number = number
        netstim.interval = interval

        return netstim


