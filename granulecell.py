# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 14:36:11 2017

@author: DanielM
"""
from genericcell import GenericCell
from neuron import h, gui

class GranuleCell(GenericCell):
    """Create a granule cell, by using the methods provided by GenericCell"""
    
    def __init__(self):
        self.all_sections = []
        # Set up the sections with topology
        self.mk_sections(2,4)
        # Set up geometry of sections
        self.mk_geometry(16.8, 16.8, [[50,150,150,150], [50,150,150,150]],[[3,3,3,3],[3,3,3,3]])
        # The Mechanisms that are equal in all sections
        for x in self.all_sections:
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
            x[0].gncabar_nca = 0.003
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
            # x.cao_ccanl = 2.0 defined in the soltesz model but not available here
        #self.pp_stim()

    def pp_stim(self):
        self.netstim = h.NetStim()
        self.netstim.interval = 100
        self.netstim.number = 1
        self.netstim.start = 1000
        self.synapse = h.Exp2Sid(self.dendrites[0][3](0.5))
        self.synapse.tau1 = 1.5
        self.synapse.tau2 = 5.5
        self.synapse.e = 0
        self.synapse.sid = 0
        #netstim.interval = 100
        #netstim.number = 1
        #netstim.start = 5
        netcon = h.NetCon(self.netstim, self.synapse, 10,1,1)
        self.netcon = netcon
        
        
        
        
        




