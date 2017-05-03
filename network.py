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
        return self.connections[name]

    def connect_cells(self, pre_type, post_type, target_sections, divergence,
                      tau1, tau2, g_max, e, thr, delay, weight,
                      name = 'default connection'):
        """Connect cells of type pre_type with cells of post_type. Return a 
        connection object that contains references to all cells, all synapses
        and all netcons involved in the connection. Right now the synapse used
        is a Exp2Sid synapse
        """
        
        conn = self.make_connection(name)

        num_synapses = len(self.cells[pre_type]) * divergence
        target_list = [x.all_sections[y] for x in self.cells[post_type] for y in target_sections]

        #Randomly select num_synapses samples out of the target list
        #rnd_sections = reservoir_sampling(target_list, num_synapses)

        #Connect a synapse to each randomly selected section
        #for x in rnd_sections:
        #    conn.synapses.append(h.Exp2Sid(x(0.5)))

        for x in range(num_synapses):
            rnd_int = random.randint(0, len(target_list) - 1)
            conn.synapses.append(h.Exp2Sid(target_list[rnd_int](0.5)))
        
        for x in conn.synapses:
            x.tau1 = tau1
            x.tau2 = tau2
            x.g = g_max
            x.e = e
        
        #Now connect a cell of pre_type to each synapse
        if (len(self.cells[pre_type]) * divergence) != len(conn.synapses):
            raise ValueError("Number of cells not equal number of synapses")
        for idx, x in enumerate(self.cells[pre_type]):
            conn.netcons.append(h.NetCon(x.soma(0.5)._ref_v, conn.synapses[idx], thr, delay, weight, sec = x.soma))
            
        self.connections[name] = conn
        
        
    def _make_synapses(self, synapse_type, cell_type, section, connection):
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

    def connect_cells_depracated(self, pre, connection, pre_divergence = 1, thr = 10, delay = 1, weight = 0):
        synapses = range(len(connection.synapses))
        if not hasattr(connection, 'netcons'):
            connection.netcons = []
        for x in self.cells[pre]:
            for y in range(pre_divergence):
                random_int = random.randint(0, len(synapses)-1)
                netcon = h.NetCon(x.soma(0.5)._ref_v, connection.synapses[synapses[random_int]], thr, delay, weight, sec = x.soma) 
                connection.netcons.append(netcon)

    def current_clamp(self, cell_type, amp, dur, delay):
        rnd_int = random.randint(0,len(self.cells[cell_type]) - 1)
        soma_v_vec, t_vec = self.cells[cell_type][rnd_int]._current_clamp_soma(amp, dur, delay)
        return soma_v_vec, t_vec

    def voltage_recording(self, cell_type):
        rnd_int = random.randint(0,len(self.cells[cell_type]) - 1)
        soma_v_vec, t_vec = self.cells[cell_type][rnd_int]._voltage_recording()
        return soma_v_vec, t_vec
        
    def run_network(self, tstop = 1000, dt = 1):
        h.tstop = tstop
        h.run()

class _Connection(object):
    """A Connection is a collection of objects that implement a type of
    connectivity. It collects the synapses, maps them to their target and their
    sources.
    """
    
    def __init__(self, name):
        self.name = name
        self.synapses = []
        self.netcons = []
        
    def add_synapses(self, synapses):
        if not hasattr(self, 'synapses'):
            self.synapses = synapses
        else:
            for x in synapses:
                self.synapses.append(x)
        
    
    
"""HELPERS"""

def reservoir_sampling(iterator, K):
    result = []
    N = 0

    for item in iterator:
        N += 1
        if len( result ) < K:
            result.append( item )
        else:
            s = int(random.random() * N)
            if s < K:
                result[ s ] = item

    return result
    
    
    
    
    
    
    
    
    
    