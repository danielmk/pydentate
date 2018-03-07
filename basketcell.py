# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:38:35 2017

@author: DanielM
"""

from ouropy.genneuron import GenNeuron
import ouropy.parameters as params

class BasketCell(GenNeuron):
    name = "BasketCell"
    def __init__(self, name = None):
        
        self.name = name
        self.all_secs = []
        self.dendrites = []

        # Make soma
        self.mk_soma(name='soma', diam=15, L=20)

        # Make dendrites
        self.mk_dendrite(4, dend_name = 'basal_1', sec_names = ['proxd', 'mid1', 'mid2', 'ddend'], diam = [4,3,2,1], L = [50,50,50,50], soma_loc = 0)
        self.mk_dendrite(4, dend_name = 'basal_2', sec_names = ['proxd', 'mid1', 'mid2', 'ddend'], diam = [4,3,2,1], L = [50,50,50,50], soma_loc = 0)
        self.mk_dendrite(4, dend_name = 'apical_1', sec_names = ['proxd', 'mid1', 'mid2', 'ddend'], diam = [4,3,2,1], L = [75,75,75,75], soma_loc = 1.0)
        self.mk_dendrite(4, dend_name = 'apical_2', sec_names = ['proxd', 'mid1', 'mid2', 'ddend'], diam = [4,3,2,1], L = [75,75,75,75], soma_loc = 1.0)

        parameters = params.read_parameters('basketcellparams.txt')

        self.insert_mechs(parameters)
        
        
"""        self.name = name
        self.all_secs = []
        self.dendrites = []
        # Set up the sections with topology
        self.mk_soma(name = 'soma', diam = 16.8, L = 16.8)
        self.mk_sections(4,4)
        # Set up geometry of sectionss
        self.mk_geometry(15, 20, [[75,75,75,75], [75,75,75,75], [50,50,50,50], [50,50,50,50]],[[4,3,2,1],[4,3,2,1], [4,3,2,1], [4,3,2,1]])
        
        for x in self.all_sections:
            x.insert('ccanl')
            x.catau_ccanl = 10
            x.caiinf_ccanl = 5 * (10 ** (-6))
            x.insert('borgka')
            x.gkabar_borgka = 0.00015
            x.insert('nca')
            x.gncabar_nca = 0.0008
            x.insert('lca')
            x.glcabar_lca = 0.005
            x.insert('gskch')
            x.gskbar_gskch = 0.000002
            x.insert('cagk')
            x.gkbar_cagk = 0.0002
            
        self.soma.insert('ichan2')
        self.soma.gnatbar_ichan2 = 0.12
        self.soma.gkfbar_ichan2 = 0.013
        self.soma.gl_ichan2 = 0.00018
        self.soma.cm = 1.4
    
        for x in self.dendrites:
        #Mechanisms of the proxd which is in [0]
            x[0].insert('ichan2')
            x[0].gnatbar_ichan2 = 0.12
            x[0].gkfbar_ichan2 = 0.013
            x[0].gl_ichan2 = 0.00018
            x[0].cm = 1.4
            
        for x in self.dendrites:
        #Mechanisms of the mid1 which is in [0]
            x[1].insert('ichan2')
            x[1].gnatbar_ichan2 = 0.12
            x[1].gkfbar_ichan2 = 0.013
            x[1].gl_ichan2 = 0.00018
            x[1].cm = 1.4
            
        for x in self.dendrites:
        #Mechanisms of the mid2 which is in [0]
            x[2].insert('ichan2')
            x[2].gnatbar_ichan2 = 0.0
            x[2].gkfbar_ichan2 = 0.0
            x[2].gl_ichan2 = 0.00018
            x[2].cm = 1.4
            
        for x in self.dendrites:
        #Mechanisms of the ddend which is in [0]
            x[3].insert('ichan2')
            x[3].gnatbar_ichan2 = 0.0
            x[3].gkfbar_ichan2 = 0.0
            x[3].gl_ichan2 = 0.00018
            x[3].cm = 1.4
            
        for x in self.all_sections:
            x.Ra = 100
            x.enat = 55
            x.ekf = -90
            x.ek = -90
            x.elca = 130
            x.esk = -90
            x.el_ichan2 = -60.06
            # x.cao_ccanl = 2.0 Defined in the soltesz model but not available here
"""