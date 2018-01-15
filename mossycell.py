# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 15:16:40 2017

@author: DanielM
"""

from ouropy.genneuron import GenNeuron
import ouropy.parameters as params

class MossyCell(GenNeuron):

    def __init__(self, name = None):
        self.name = name
        self.all_secs = []
        self.dendrites = []
        
        self.mk_soma(name = 'soma', diam = 20, L = 20)
        # Set up the sections with topology
        #NOTE! There is a discrepancy here between the diameter of the third section of the dendrite. Publication: 25; Code: 2.5 The number in the code seems to be more physiological
        self.mk_dendrite(4, dend_name = 'dend_1', sec_names = ['proxd', 'mid1d', 'mid2d', 'dd'], diam = [5.78,4,2.5,1], L = [50,50,50,50], soma_loc = 1)
        self.mk_dendrite(4, dend_name = 'dend_2', sec_names = ['proxd', 'mid1d', 'mid2d', 'dd'], diam = [5.78,4,2.5,1], L = [50,50,50,50], soma_loc = 1)
        self.mk_dendrite(4, dend_name = 'dend_3', sec_names = ['proxd', 'mid1d', 'mid2d', 'dd'], diam = [5.78,4,2.5,1], L = [50,50,50,50], soma_loc = 0)
        self.mk_dendrite(4, dend_name = 'dend_4', sec_names = ['proxd', 'mid1d', 'mid2d', 'dd'], diam = [5.78,4,2.5,1], L = [50,50,50,50], soma_loc = 0)

        parameters = params.read_parameters('mossycellparams_cat.txt')

        self.insert_mechs(parameters)
































"""
        for x in self.all_sections:
            x.insert('ccanl')
            x.catau_ccanl = 10
            x.caiinf_ccanl = 5 * (10 ** (-6))
            x.insert('borgka')
            x.gkabar_borgka = 0.00001
            x.insert('nca')
            x.gncabar_nca = 0.00008
            x.insert('lca')
            x.glcabar_lca = 0.0006
            x.insert('gskch')
            x.gskbar_gskch = 0.016
            x.insert('cagk')
            x.gkbar_cagk = 0.0165
            x.insert('hyperde3')
            x.ghyfbar_hyperde3 = 0.000005
            x.ghysbar_hyperde3 = 0.000005

        self.soma.insert('ichan2')
        self.soma.gnatbar_ichan2 = 0.12
        self.soma.gkfbar_ichan2 = 0.0005
        self.soma.gl_ichan2 = 0.000011
        self.soma.cm = 0.6

        for x in self.dendrites:
            x[0].insert('ichan2')
            x[0].gnatbar_ichan2 = 0.12
            x[0].gkfbar_ichan2 = 0.0005
            x[0].gl_ichan2 = 0.000044
            x[0].cm = 2.4

        for x in self.dendrites:
            for y in x[1:]:
                y.insert('ichan2')
                y.gnatbar_ichan2 = 0.0
                y.gkfbar_ichan2 = 0.00
                y.gl_ichan2 = 0.000044
                y.cm = 2.4
                
        for x in self.all_sections:
            x.Ra = 100
            x.enat = 55
            x.ekf = -90
            x.ek = -90
            x.esk = -90
            x.elca = 130
            x.ehyf = -40
            x.ehys = -40
            x.el_ichan2 = -59
            #x.cao_ccanl = 2 Defined in the soltesz model but not available here
                """