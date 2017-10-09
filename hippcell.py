# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 10:42:49 2017

@author: DanielM
"""

from ouropy.genneuron import GenNeuron
import ouropy.parameters as params

class HippCell(GenNeuron):
    def __init__(self, name = None):
        self.name = name
        self.all_secs = []
        self.dendrites = []
        
        self.mk_soma(name='soma', diam=10, L=20)
        
        self.mk_dendrite(3, dend_name = 'short_1', sec_names = ['proxd', 'midd', 'distd'], diam = [3,2,1], L = [50,50,50], soma_loc = 0.0)
        self.mk_dendrite(3, dend_name = 'short_2', sec_names = ['proxd', 'midd', 'distd'], diam = [3,2,1], L = [50,50,50], soma_loc = 0.0)
        self.mk_dendrite(3, dend_name = 'long_1', sec_names = ['proxd', 'midd', 'distd'], diam = [3,2,1], L = [75,75,75], soma_loc = 1.0)
        self.mk_dendrite(3, dend_name = 'long_2', sec_names = ['proxd', 'midd', 'distd'], diam = [3,2,1], L = [75,75,75], soma_loc = 1.0)

        parameters = params.read_parameters('hippcellparams.txt')

        self.parameters = parameters

        self.insert_mechs(parameters)

"""self.all_sections = []
            # Set up the sections with topology
            self.mk_sections(4,3)
            # Set up geometry of sectionss
            self.mk_geometry(10, 20, [[75,75,75], [75,75,75], [50,50,50], [50,50,50]],[[3,2,1],[3,2,1], [3,2,1], [3,2,1]])

            for x in self.all_sections:
                x.insert('ccanl')
                x.catau_ccanl = 10
                x.caiinf_ccanl = 5 * (10 ** (-6))
                x.insert('borgka')
                x.gkabar_borgka = 0.0008
                x.insert('nca')
                x.gncabar_nca = 0.0
                x.insert('lca')
                x.glcabar_lca = 0.0015
                x.insert('gskch')
                x.gskbar_gskch = 0.003
                x.insert('cagk')
                x.gkbar_cagk = 0.003
                x.insert('hyperde3')
                x.ghyfbar_hyperde3 = 0.000015
                x.ghysbar_hyperde3 = 0.000015
                
            self.soma.insert('ichan2')
            self.soma.gnatbar_ichan2 = 0.2
            self.soma.gkfbar_ichan2 = 0.006
            self.soma.gl_ichan2 = 0.000036
            self.soma.cm = 1.1
        
            for x in self.dendrites:
            #Mechanisms of the proxd which is in [0]
                x[0].insert('ichan2')
                x[0].gnatbar_ichan2 = 0.2
                x[0].gkfbar_ichan2 = 0.006
                x[0].gl_ichan2 = 0.000036
                x[0].cm = 1.1
                
            for x in self.dendrites:
            #Mechanisms of the midd, distd which is in [1] and [2]
                for y in x[1:]:
                    y.insert('ichan2')
                    y.gnatbar_ichan2 = 0.0
                    y.gkfbar_ichan2 = 0.0
                    y.gl_ichan2 = 0.000036
                    y.cm = 1.1
                
            for x in self.all_sections:
                x.Ra = 100
                x.enat = 55
                x.ekf = -90
                x.ek = -90
                x.elca = 130
                x.esk = -90
                x.el_ichan2 = -70.45
                x.ehyf = -40
                x.ehys = -40
                # x.cao_ccanl = 2.0 defined in the soltesz model but not available here"""