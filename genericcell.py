# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:33:43 2017

@author: DanielM
"""

from neuron import h, gui

class GenericCell(object):
    """This is a generic cell which implements the five compontents of neuron:
    Sections - the cell sections
    Topology - the connectivity of the sections
    Geometry - The 3D location of the sections
    Biophysics - The ionic channels and membrane properties of the sections
    Synapses - Optional list of synapses onto the cell.
    It implements those in a generic way so a class of a specific cell type
    should not need to overwrite the components, but could just call the
    methods to initialize it's parameters.
    """

    def __init__(self):
        """Generic __init__ mostly for debugging"""
    def mk_sections(self, num_dend, num_seg):
        """Sets up the Sections WITH Topology"""
        self.soma = h.Section(name='soma')
        self.all_sections.append(self.soma)
        self.dendrites = []
        self.num_dend = num_dend
        self.i = 0  # Counter for __iter__

        if hasattr(num_seg, '__iter'):
            for curr_dend in range(num_dend):
                for curr_num_segs in num_seg:
                    self.dendrites.append(Dendrite(n_segs = curr_num_segs))
        else:
            for curr_dend in range(num_dend):
                self.dendrites.append(Dendrite(n_segs = num_seg))
        for x in self.dendrites:
            for y in x:
                self.all_sections.append(y)

        for x in self.dendrites:
            x.connect_segments(self.soma)

    def mk_geometry(self, soma_diam, soma_L, dend_L, dend_diam):
        """Sets up the geometry of the sections
        dend_L - scalar or iterable matching the number of dendrites
        dend_diam - scalar or iterable matching the number of dendrites
        """
        self.soma.diam = soma_diam
        self.soma.L = soma_L

        if len(dend_L) != self.num_dend:
            raise ValueError("dend_L must match the number of dendrites or be a scalar")
        if len(dend_diam) != self.num_dend:
            raise ValueError("dend_diam must match the number of dendrites or be a scalar")
        for dend_idx, x in enumerate(self.dendrites):
            if len(self.dendrites[dend_idx]) != len(dend_L[dend_idx]):
                raise ValueError
            if len(self.dendrites[dend_idx]) != len(dend_diam[dend_idx]):
                raise ValueError
            for seg_idx, y in enumerate(x):
                y.L = dend_L[dend_idx][seg_idx]
                y.diam = dend_diam[dend_idx][seg_idx]
                
    def connect_post(self, source, synapse, thr = 10, delay = 1, weight = 0):
        if not hasattr(self, 'synapses'):
            self.synapses = []
        if not hasattr(self, 'netcons'):
            self.netcons = []

        self.synapses.append(synapse)
        netcon = h.NetCon(source, synapse, thr, delay, weight)
        self.netcons.append(netcon)

        return netcon

    def somatic_recording(self):
        self.stim = h.IClamp(self.soma(1))
        self.stim.delay = 500
        self.stim.dur = 500
        self.stim.amp = 300
        
        soma_v_vec = h.Vector()
        t_vec = h.Vector()
        soma_v_vec.record(self.soma(1)._ref_v)
        t_vec.record(h._ref_t)

        return soma_v_vec, t_vec

    def __iter__(self):
        return self

    def next(self):
        if self.i < (len(self.segs)):
            i = self.i
            self.i += 1
            return self.all_sections[i]
        else:
            self.i = 0
            raise StopIteration()
            

    

class Dendrite(object):
    def __init__(self, name = "dendrite", n_segs= 0, diam = None, L = None):
        self.name = name
        self.segs = []
        self.i = 0
        self.mk_segments(n_segs = n_segs)
        if bool(diam):
            self.set_diam(diam)
        if bool(L):
            self.set_L(L)

    def mk_segments(self, n_segs = 1):
        self.segs = []
        for curr_n in range(n_segs):
            self.segs.append(h.Section(name = self.name + 'seg' + str(curr_n)))

    def connect_segments(self, soma):
        for idx, curr_seg in enumerate(self.segs):
            if idx == 0:
                curr_seg.connect(soma(1))
            else:
                curr_seg.connect(self.segs[idx -1](1))

    def set_diam(self, diam):
        if not bool(self.segs):
            raise Warning("Can't set diameter before segments are made")
            return
        if hasattr(diam, '__iter__'):
            if len(diam) != len(self.segs):
                raise Warning("List of diameters does not fit number of segments")
                return
            for idx, curr_seg in enumerate(self.segs):
                curr_seg.diam = diam[idx]
            return
        else:
            for curr_seg in self.segs:
                curr_seg.diam = diam

    def set_L(self, L):
        if not bool(self.segs):
            raise Warning("Can't set L before segments are made")
            return
        if hasattr(L, '__iter__'):
            if len(L) != len(self.segs):
                raise Warning("List of diameters does not fit number of segments")
                return
            for idx, curr_seg in enumerate(self.segs):
                curr_seg.L = L[idx]
            return
        else:
            for curr_seg in self.segs:
                curr_seg.L = L

    def __iter__(self):
        return self

    def next(self):
        if self.i < (len(self.segs)):
            i = self.i
            self.i += 1
            return self.segs[i]
        else:
            self.i = 0
            raise StopIteration()

    def __getitem__(self, key):
        return self.segs[key]

    def __len__(self):
        return len(self.segs)

    
    
    