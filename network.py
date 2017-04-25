# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 11:16:23 2017

@author: DanielM
"""

from neuron import h, gui
from granulecell import GranuleCell
from mossycell import MossyCell
from basketcell import BasketCell
from hippcell import HippCell
import random
import copy

class Network(object):
    """The Network object can create a number of cells, connect them with
    artificial sources or connect the cells with each other. It also keeps
    track of the different cells, synapses and connections and other properties
    of the network. In general: any property that is not an intrinsic property
    of a cell is a property of the network and thus managed here.
    """
    def __init__(self, celltypes, cellnums):
        self.cells = {}
        self.connections = {}
        for idx, x in enumerate(celltypes):
            self.cells[x] = []
            for y in range(cellnums[idx]):
                self.cells[x].append(x())

    def make_connection(self, name):
        self.connections[name] = _Connection(name)        

    def make_synapses(self, synapse_type, cell_type, section, connection):
        synapse_list_tmp = []
        for x in self.cells[cell_type]:
            if hasattr(section, '__iter__'):
                for y in section:
                    curr_synapse = synapse_type(x.all_sections[y](0.5))
                    synapse_list_tmp.append(curr_synapse)
            else:
                curr_synapse = synapse_type(x.all_sections[section](0.5))
                synapse_list_tmp.append(curr_synapse)
        
        connection.add_synapses(synapse_list_tmp)

    def connect_synapses_artificial(self, source, connection, thr = 10, delay = 1, weight = 0):
        """Connect synapses with an artificial source_type"""
        if not hasattr(connection, 'sources'):
            connection.sources = []
            connection.sources.append(source)
        else:
            connection.sources.append(source)
        
        if not hasattr(connection, 'netcons'):
            connection.netcons = []
        
        for x in connection.synapses:
            netcon = h.NetCon(source, x, thr, delay, weight)
            connection.netcons.append(netcon)

    def connect_cells(self, pre, connection, pre_divergence = 1, thr = 10, delay = 1, weight = 0):
        synapses = range(len(connection.synapses))
        if not hasattr(connection, 'netcons'):
            connection.netcons = []
        for x in self.cells[pre]:
            for y in range(pre_divergence):
                random_int = random.randint(0, len(synapses)-1)
                netcon = h.NetCon(x.soma(0.5)._ref_v, connection.synapses[synapses[random_int]], thr, delay, weight, sec = x.soma) 
                connection.netcons.append(netcon)
                del synapses[random_int]


    def test_recording(self, cell_type):

        rnd_int = random.randint(0,len(self.cells[cell_type]) - 1)

        soma_v_vec, t_vec = self.cells[cell_type][rnd_int].somatic_recording()

        return soma_v_vec, t_vec
        
    def run_network(self):
        h.tstop = 2000
        print(h.tstop)
        h.run()

class _Connection(object):
    """A Connection is a collection of objects that implement a type of
    connectivity. It collects the synapses, maps them to their target and their
    sources.
    """
    
    def __init__(self, name):
        self.name = name
        
    def add_synapses(self, synapses):
        if not hasattr(self, 'synapses'):
            self.synapses = synapses
        else:
            for x in synapses:
                self.synapses.append(x)
        
    
    
    
    
    
    
    
    
    
    
    
    
    