# -*- coding: utf-8 -*-
"""
Created on Wed Apr 05 12:36:11 2017

@author: DanielM
"""

"""
This is a fully wired network that functions with 50000 GCs and 1 PP input
Auto init and run
"""

import numpy as np
from neuron import h, gui
import matplotlib.pyplot as plt

h.nrn_load_dll("C:\\nrn\\dentate_gyrus_python_translate\\nrnmech.dll")

# Set up some simulation parameters

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
        self.stim.amp = .5
    
        soma_v_vec = h.Vector()
        t_vec = h.Vector()
        soma_v_vec.record(self.soma(0.5)._ref_v)
        t_vec.record(h._ref_t)
        
        print("Somatic recording set up")
        return soma_v_vec, t_vec
        
    def simulate(self):
        print("Starting simulation")
        h.tstop = 1500
        h.run()
        print("Stimulation finished")
        

class GranuleCell(GenericCell):
    """Create a granule cell, by using the methods provided by GenericCell"""
    
    def __init__(self):
        self.all_sections = []
        # Set up the sections with topology
        self.mk_sections(2,4)
        # Set up geometry of sectionss
        self.mk_geometry(16.8, 16.8, [[50,150,150,150], [50,150,150,150]],[[3,3,3,3],[3,3,3,3]])
        # The Mechanisms that are equal in all sections
        for x in self.all_sections:
            print(x.name())
            x.insert('ccanl')
            x.catau_ccanl = 10
            x.caiinf_ccanl = 5 * (10 ** (-6))
            x.Ra = 210

            
        # The Mechanisms in the soma
        self.soma.insert('ichan2')
        self.soma.gnatbar_ichan2 = 0.12
        self.soma.gkfbar_ichan2  =0.016
        self.soma.gksbar_ichan2 = 0.006
        self.soma.insert('borgka')
        self.soma.gkabar_borgka = 0.012
        self.soma.insert('nca')
        self.soma.gncabar_nca = 0.002
        self.soma.insert('lca')
        self.soma.glcabar_lca = 0.005
        self.soma.insert('cat')
        self.soma.gcatbar_cat = 0.000037
        self.soma.insert('gskch')
        self.soma.gskbar_gskch = 0.001
        self.soma.insert('cagk')
        self.soma.gkbar_cagk = 0.0006
        self.soma.gl_ichan2 = 0.00004
        self.soma.cm = 1
        
        #The Dendrite mechanisms
        for x in self.dendrites:
            #Mechanisms of the gcldend which is in [0]
            x[0].insert('ichan2')
            x[0].gnatbar_ichan2 = 0.018
            x[0].gkfbar_ichan2 = 0.004
            x[0].gksbar_ichan2 = 0.006
            x[0].insert('nca')
            x[0].gncabar_nca = 0.004
            x[0].insert('lca')
            x[0].glcabar_lca = 0.0075
            x[0].insert('cat')
            x[0].gcatbar_cat = 0.000075
            x[0].insert('gskch')
            x[0].gskbar_gskch = 0.0004
            x[0].insert('cagk')
            x[0].gkbar_cagk = 0.0006
            x[0].gl_ichan2 = 0.00004
            x[0].cm = 1
        
        for x in self.dendrites:
            #Mechanisms of the proximal dend which is in [1]
            x[1].insert('ichan2')
            x[1].gnatbar_ichan2 = 0.013
            x[1].gkfbar_ichan2 = 0.004
            x[1].gksbar_ichan2 = 0.006
            x[1].insert('nca')
            x[1].gncabar_nca = 0.001
            x[1].insert('lca')
            x[1].glcabar_lca = 0.0075
            x[1].insert('cat')
            x[1].gcatbar_cat = 0.00025
            x[1].insert('gskch')
            x[1].gskbar_gskch = 0.0002
            x[1].insert('cagk')
            x[1].gkbar_cagk = 0.001
            x[1].gl_ichan2 = 0.000063
            x[1].cm = 1.6
            
        for x in self.dendrites:
            #Mechanisms of the medial dend which is in [2]
            x[2].insert('ichan2')
            x[2].gnatbar_ichan2 = 0.008
            x[2].gkfbar_ichan2 = 0.001
            x[2].gksbar_ichan2 = 0.006
            x[2].insert('nca')
            x[2].gncabar_nca = 0.001
            x[2].insert('lca')
            x[2].glcabar_lca = 0.0005
            x[2].insert('cat')
            x[2].gcatbar_cat = 0.0005
            x[2].insert('gskch')
            x[2].gskbar_gskch =0.0
            x[2].insert('cagk')
            x[2].gkbar_cagk = 0.0024
            x[2].gl_ichan2 = 0.000063
            x[2].cm = 1.6
            
        for x in self.dendrites:
            #Mechanisms of the distal dend which is in [3]
            x[3].insert('ichan2')
            x[3].gnatbar_ichan2 = 0.0
            x[3].gkfbar_ichan2 = 0.001
            x[3].gksbar_ichan2 = 0.008
            x[3].insert('nca')
            x[3].gncabar_nca = 0.001
            x[3].insert('lca')
            x[3].glcabar_lca = 0.0
            x[3].insert('cat')
            x[3].gcatbar_cat = 0.001
            x[3].insert('gskch')
            x[3].gskbar_gskch = 0.0
            x[3].insert('cagk')
            x[3].gkbar_cagk = 0.0024
            x[3].gl_ichan2 = 0.000063
            x[3].cm = 1.6
        
        for x in self.all_sections:
            x.enat = 45
            x.ekf = -90
            x.eks = -90
            x.ek = -90
            x.elca = 130
            x.etca = 130
            x.esk = -90
            x.el_ichan2 = -70
            x.cai_ccanl = 2.0
            
        print("GranuleCell initialized")




