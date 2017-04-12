# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:33:43 2017

@author: DanielM
"""

from neuron import h, gui

h.nrn_load_dll("C:\\nrn\\dentate_gyrus_python_translate\\nrnmech_new.dll")

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
        self.mk_sections(2,4)
        self.mk_geometry(16.8, 16.8, [[50,150,150,150], [50,150,150,150]],[[3,3,3,3],[3,3,3,3]])
        self.all_sections = []
    def mk_sections(self, num_dend, num_sec):
        """Sets up the Sections WITH Topology"""
        self.soma = h.Section(name = 'soma')
        self.all_sections.append(self.soma)
        self.dendrites = []
        
        for x in range(num_dend):
            self.dendrites.append([])
            for y in range(num_sec):
                curr_seg = h.Section(name = "dend_" + str(x) + "_sec_" + str(y))
                if y == 0:
                    curr_seg.connect(self.soma)
                else:
                    curr_seg.connect(prev_seg)
                prev_seg = curr_seg
                self.dendrites[x].append(curr_seg)
                self.all_sections.append(curr_seg)
            self.num_dend = len(self.dendrites)
    #def mk_topology():
                
    def mk_geometry(self, soma_diam,soma_L, dend_L, dend_diam):
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

    def somatic_recording(self):
        self.stim = h.IClamp(self.soma(1))
        self.stim.delay = 500
        self.stim.dur = 500
        self.stim.amp = .3
    
        soma_v_vec = h.Vector()
        t_vec = h.Vector()
        soma_v_vec.record(self.soma(0.5)._ref_v)
        t_vec.record(h._ref_t)
        
        print("Somatic recording set up")
        return soma_v_vec, t_vec
        
    def simulate(self):
        print("Starting simulation")
        h.tstop = 10000
        h.run()
        print("Stimulation finished")