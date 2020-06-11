# -*- coding: utf-8 -*-
"""
This module implements the GenNetwork class, which implements generic network
logic. Also implements Population and GenConnection classes

@author: DanielM
"""

from neuron import h
import random
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import os
import shelve
import scipy.stats as stats


class GenConnection(object):
    def __init__(self):
        pass

    def get_description(self):
        """Return a descriptive string for the connection"""
        name = self.pre_pop.name + ' to ' + self.post_pop.name + '\n'
        pre_cell_targets = '\n'.join([str(x) for x in self.pre_cell_targets])
        return name + pre_cell_targets

    def get_name(self):
        if type(self.pre_pop) == str:
            return self.pre_pop + ' to ' + str(self.post_pop)
        else:
            return str(self.pre_pop) + ' to ' + str(self.post_pop)

    def get_properties(self):
        """Get the and make them suitable for pickling"""
        properties = {'name': self.get_name(),
                      'init_parameters': self.init_parameters,
                      'pre_cell_targets': self.pre_cell_targets}
        properties['init_parameters']['post_pop'] = str(properties['init_parameters']['post_pop'])
        properties['init_parameters']['self'] = str(properties['init_parameters']['self'])
        try:
             properties['init_parameters']['pre_pop'] = str(properties['init_parameters']['pre_pop'])
        except:
            pass
        return {self.get_name(): properties}


class tmgsynConnection(GenConnection):

    def __init__(self, pre_pop, post_pop,
                 target_pool, target_segs, divergence,
                 tau_1, tau_facil, U, tau_rec, e, thr, delay, weight):
        """Create a connection with tmgsyn as published by Tsodyks, Pawelzik &
        Markram, 1998.
        The tmgsyn is a dynamic three state implicit resource synapse model.
        The response onset is instantaneous and the decay is exponential.
        It combines a frequency dependent depression and a facilitation
        mechanism that both depend on independent time constants.
        The synaptic targets are chosen by proximity, that is the target pool
        are the cells in closest proximity.

        Parameters
        ----------
        pre_pop - gennetwork.Population
            The presynaptic population
        post_pop - gennetwork.Population
            the postsynaptic population
        target_pool - int
            the number of cells in the target pool
        target_segs - str
            the name of the segments that are possible synaptic targets at the
            postsynaptic population
        divergence - int
            divergence in absolute terms, that is the number of synapses each
            presynaptic cell forms
        tau_1 - numeric
            the time constant of synaptic decay. conforms to the transition
            from the active to the inactive resource state. units of time as in
            neuron standard units
        tau_facil - numeric
            the time constant of facilitation decay. this essentially creates
            the frequency dependence. set to 0 for no facilitation.
        U - numeric
            maximum of postsynaptic response. has to be considered together
            with the weight from the netcon.
        tau_rec - numeric
            time constant of recovery from inactive for recovered state.
            gives depression since inactive resources do not contribute to
            postsynaptic signal. set to 0 for no depression.
        e - numeric
            equilibrium potential of the postsynaptic conductance
        thr - numeric
            threshold for synaptic event at the presynaptic source
        delay - numeric
            delay between presynaptic signal and onset of postsynaptic signal
        weight - numeric
            weight for the netcon object connecting source and target

        Returns
        -------
        None

        Use Cases
        ---------
        >>> tmgsynConnection(nw.population[0], nw.population[1],
                             3, 'prox', 1, 6.0, 0, 0.04, 0, 0, 10, 3, 0)
        A non-facilitating, non-depressing excitatory connection.

        """
        self.init_parameters = locals()
        self.pre_pop = pre_pop
        self.post_pop = post_pop
        pre_pop.add_connection(self)
        post_pop.add_connection(self)
        pre_pop_rad = (np.arange(pre_pop.get_cell_number(), dtype=float) /
                       pre_pop.get_cell_number()) * (2*np.pi)
        post_pop_rad = (np.arange(post_pop.get_cell_number(), dtype=float) /
                        post_pop.get_cell_number()) * (2*np.pi)

        pre_pop_pos = pos(pre_pop_rad)
        post_pop_pos = pos(post_pop_rad)
        pre_cell_target = []
        synapses = []
        netcons = []
        conductances = []

        for idx, curr_cell_pos in enumerate(pre_pop_pos):

            curr_dist = []
            for post_cell_pos in post_pop_pos:
                curr_dist.append(euclidian_dist(curr_cell_pos, post_cell_pos))

            sort_idc = np.argsort(curr_dist)
            closest_cells = sort_idc[0:target_pool]
            picked_cells = np.random.choice(closest_cells,
                                            divergence,
                                            replace=False)
            pre_cell_target.append(picked_cells)
            for tar_c in picked_cells:

                curr_syns = []
                curr_netcons = []
                curr_conductances = []

                curr_seg_pool = post_pop[tar_c].get_segs_by_name(target_segs)
                chosen_seg = np.random.choice(curr_seg_pool)
                for seg in chosen_seg:
                    curr_syn = h.tmgsyn(chosen_seg(0.5))
                    curr_syn.tau_1 = tau_1
                    curr_syn.tau_facil = tau_facil
                    curr_syn.U = U
                    curr_syn.e = e
                    curr_syn.tau_rec = tau_rec
                    curr_syns.append(curr_syn)
                    curr_netcon = h.NetCon(pre_pop[idx].soma(0.5)._ref_v,
                                           curr_syn, thr, delay,
                                           weight, sec=pre_pop[idx].soma)
                    curr_gvec = h.Vector()
                    curr_gvec.record(curr_syn._ref_g)
                    curr_conductances.append(curr_gvec)
                    curr_netcons.append(curr_netcon)
                    netcons.append(curr_netcons)
                    synapses.append(curr_syns)
            conductances.append(curr_conductances)
        self.conductances = conductances
        self.netcons = netcons
        self.pre_cell_targets = np.array(pre_cell_target)
        self.synapses = synapses


class pyr2pyrConnCA3(GenConnection):

    def __init__(self, pre_pop, post_pop,
                 target_pool, target_segs, divergence, thr, weight,
                 initW, Max, Min, Delay, Ratio, ConPattern, Lambda1, Lambda2,
                 threshold1, threshold2, tauD1, d1, tauD2, d2, tauF, f, erev,
                 bACH, aDA, bDA, wACH, Alpha_ampa, Beta_ampa, Cdur_ampa,
                 gbar_ampa, Alpha_nmda, Beta_nmda, Cdur_nmda, gbar_nmda,
                 GaussA):
        """TODO DOC

        Parameters
        ----------
        pre_pop - gennetwork.Population
            The presynaptic population
        post_pop - gennetwork.Population
            the postsynaptic population
        target_pool - int
            the number of cells in the target pool
        target_segs - str
            the name of the segments that are possible synaptic targets at the
            postsynaptic population
        divergence - int
            divergence in absolute terms, that is the number of synapses each
            presynaptic cell forms
        tau_1 - numeric
            the time constant of synaptic decay. conforms to the transition
            from the active to the inactive resource state. units of time as in
            neuron standard units
        tau_facil - numeric
            the time constant of facilitation decay. this essentially creates
            the frequency dependence. set to 0 for no facilitation.
        U - numeric
            maximum of postsynaptic response. has to be considered together
            with the weight from the netcon.
        tau_rec - numeric
            time constant of recovery from inactive for recovered state.
            gives depression since inactive resources do not contribute to
            postsynaptic signal. set to 0 for no depression.
        e - numeric
            equilibrium potential of the postsynaptic conductance
        thr - numeric
            threshold for synaptic event at the presynaptic source
        delay - numeric
            delay between presynaptic signal and onset of postsynaptic signal
        weight - numeric
            weight for the netcon object connecting source and target

        Returns
        -------
        None

        Use Cases
        ---------
        >>> tmgsynConnection(nw.population[0], nw.population[1],
                             3, 'prox', 1, 6.0, 0, 0.04, 0, 0, 10, 3, 0)
        A non-facilitating, non-depressing excitatory connection.

        """
        self.init_parameters = locals()
        self.pre_pop = pre_pop
        self.post_pop = post_pop
        pre_pop.add_connection(self)
        post_pop.add_connection(self)
        pre_pop_rad = (np.arange(pre_pop.get_cell_number(), dtype=float) /
                       pre_pop.get_cell_number()) * (2*np.pi)
        post_pop_rad = (np.arange(post_pop.get_cell_number(), dtype=float) /
                        post_pop.get_cell_number()) * (2*np.pi)

        pre_pop_pos = pos(pre_pop_rad)
        post_pop_pos = pos(post_pop_rad)
        pre_cell_target = []
        synapses = []
        netcons = []
        conductances = []

        for idx, curr_cell_pos in enumerate(pre_pop_pos):

            curr_dist = []
            for post_cell_pos in post_pop_pos:
                curr_dist.append(euclidian_dist(curr_cell_pos, post_cell_pos))

            sort_idc = np.argsort(curr_dist)
            closest_cells = sort_idc[0:target_pool]
            picked_cells = np.random.choice(closest_cells,
                                            divergence,
                                            replace=False)
            pre_cell_target.append(picked_cells)
            for tar_c in picked_cells:

                curr_syns = []
                curr_netcons = []
                curr_conductances = []

                curr_syn = h.pyr2pyr(post_pop[tar_c].soma(0.5))
                curr_syn.initW = initW
                curr_syn.Wmax = Max
                curr_syn.Wmin = Min
                curr_syn.lambda1 = Lambda1
                curr_syn.lambda2 = Lambda2
                curr_syn.threshold1 = threshold1
                curr_syn.threshold2 = threshold2
                curr_syn.tauD1 = tauD1
                curr_syn.d1 = d1
                curr_syn.tauD2 = tauD2
                curr_syn.d2 = d2
                curr_syn.tauF = tauF
                curr_syn.f = f
                curr_syn.Erev_ampa = erev
                curr_syn.Erev_nmda = erev
                curr_syn.bACH = bACH
                curr_syn.aDA = aDA
                curr_syn.bDA = bDA
                curr_syn.wACH = wACH
                curr_syn.AlphaTmax_ampa = Alpha_ampa
                curr_syn.Beta_ampa = Beta_ampa
                curr_syn.Cdur_ampa = Cdur_ampa
                curr_syn.gbar_ampa = gbar_ampa
                curr_syn.AlphaTmax_nmda = Alpha_nmda
                curr_syn.Beta_nmda = Beta_nmda
                curr_syn.Cdur_nmda = Cdur_nmda
                curr_syn.gbar_nmda = gbar_nmda

                curr_syns.append(curr_syn)
                curr_netcon = h.NetCon(pre_pop[idx].soma(0.5)._ref_v,
                                       curr_syn, thr, Delay,
                                       weight, sec=pre_pop[idx].soma)
                #curr_gvec = h.Vector()
                #curr_gvec.record(curr_syn._ref_g)
                #curr_conductances.append(curr_gvec)
                curr_netcons.append(curr_netcon)
                netcons.append(curr_netcons)
                synapses.append(curr_syns)
            conductances.append(curr_conductances)
        self.conductances = conductances
        self.netcons = netcons
        self.pre_cell_targets = np.array(pre_cell_target)
        self.synapses = synapses
        
class inter2pyrConnCA3(GenConnection):

    def __init__(self, pre_pop, post_pop,
                 target_pool, target_segs, divergence, thr, weight,
                 initW, Max, Min, Delay, Ratio, ConPattern, Lambda1, Lambda2,
                 threshold1, threshold2, tauD1, d1, tauD2, d2, tauF, f, erev,
                 bACH, aDA, bDA, wACH, Alpha_gabaa, Beta_gabaa, Cdur_gabaa,
                 gbar_gabaa, Alpha_gabab, Beta_gabab, Cdur_gabab, gbar_gabab,
                 GaussA):
        """TODO DOC

        Parameters
        ----------
        pre_pop - gennetwork.Population
            The presynaptic population
        post_pop - gennetwork.Population
            the postsynaptic population
        target_pool - int
            the number of cells in the target pool
        target_segs - str
            the name of the segments that are possible synaptic targets at the
            postsynaptic population
        divergence - int
            divergence in absolute terms, that is the number of synapses each
            presynaptic cell forms
        tau_1 - numeric
            the time constant of synaptic decay. conforms to the transition
            from the active to the inactive resource state. units of time as in
            neuron standard units
        tau_facil - numeric
            the time constant of facilitation decay. this essentially creates
            the frequency dependence. set to 0 for no facilitation.
        U - numeric
            maximum of postsynaptic response. has to be considered together
            with the weight from the netcon.
        tau_rec - numeric
            time constant of recovery from inactive for recovered state.
            gives depression since inactive resources do not contribute to
            postsynaptic signal. set to 0 for no depression.
        e - numeric
            equilibrium potential of the postsynaptic conductance
        thr - numeric
            threshold for synaptic event at the presynaptic source
        delay - numeric
            delay between presynaptic signal and onset of postsynaptic signal
        weight - numeric
            weight for the netcon object connecting source and target

        Returns
        -------
        None

        Use Cases
        ---------
        >>> tmgsynConnection(nw.population[0], nw.population[1],
                             3, 'prox', 1, 6.0, 0, 0.04, 0, 0, 10, 3, 0)
        A non-facilitating, non-depressing excitatory connection.

        """
        self.init_parameters = locals()
        self.pre_pop = pre_pop
        self.post_pop = post_pop
        pre_pop.add_connection(self)
        post_pop.add_connection(self)
        pre_pop_rad = (np.arange(pre_pop.get_cell_number(), dtype=float) /
                       pre_pop.get_cell_number()) * (2*np.pi)
        post_pop_rad = (np.arange(post_pop.get_cell_number(), dtype=float) /
                        post_pop.get_cell_number()) * (2*np.pi)

        pre_pop_pos = pos(pre_pop_rad)
        post_pop_pos = pos(post_pop_rad)
        pre_cell_target = []
        synapses = []
        netcons = []
        conductances = []

        for idx, curr_cell_pos in enumerate(pre_pop_pos):

            curr_dist = []
            for post_cell_pos in post_pop_pos:
                curr_dist.append(euclidian_dist(curr_cell_pos, post_cell_pos))

            sort_idc = np.argsort(curr_dist)
            closest_cells = sort_idc[0:target_pool]
            picked_cells = np.random.choice(closest_cells,
                                            divergence,
                                            replace=False)
            pre_cell_target.append(picked_cells)
            for tar_c in picked_cells:

                curr_syns = []
                curr_netcons = []
                curr_conductances = []

                curr_syn = h.inter2pyr(post_pop[tar_c].soma(0.5))
                curr_syn.initW = initW
                curr_syn.Wmax = Max
                curr_syn.Wmin = Min
                curr_syn.lambda1 = Lambda1
                curr_syn.lambda2 = Lambda2
                curr_syn.threshold1 = threshold1
                curr_syn.threshold2 = threshold2
                curr_syn.tauD1 = tauD1
                curr_syn.d1 = d1
                curr_syn.tauD2 = tauD2
                curr_syn.d2 = d2
                curr_syn.tauF = tauF
                curr_syn.f = f
                curr_syn.Erev_gabaa = erev
                curr_syn.Erev_gabab = erev
                curr_syn.bACH = bACH
                curr_syn.aDA = aDA
                curr_syn.bDA = bDA
                curr_syn.wACH = wACH
                curr_syn.AlphaTmax_gabaa = Alpha_gabaa
                curr_syn.Beta_gabaa = Beta_gabaa
                curr_syn.Cdur_gabaa = Cdur_gabaa
                curr_syn.gbar_gabaa = gbar_gabaa
                curr_syn.AlphaTmax_gabab = Alpha_gabab
                curr_syn.Beta_gabab = Beta_gabab
                curr_syn.Cdur_gabab = Cdur_gabab
                curr_syn.gbar_gabab = gbar_gabab

                curr_syns.append(curr_syn)
                curr_netcon = h.NetCon(pre_pop[idx].soma(0.5)._ref_v,
                                       curr_syn, thr, Delay,
                                       weight, sec=pre_pop[idx].soma)
                #curr_gvec = h.Vector()
                #curr_gvec.record(curr_syn._ref_g)
                #curr_conductances.append(curr_gvec)
                curr_netcons.append(curr_netcon)
                netcons.append(curr_netcons)
                synapses.append(curr_syns)
            conductances.append(curr_conductances)
        self.conductances = conductances
        self.netcons = netcons
        self.pre_cell_targets = np.array(pre_cell_target)
        self.synapses = synapses


class tmgsynConnectionExponentialProb(GenConnection):

    def __init__(self, pre_pop, post_pop,
                 scale, target_segs, divergence,
                 tau_1, tau_facil, U, tau_rec, e, thr, delay, weight):
        """Create a connection with tmgsyn as published by Tsodyks, Pawelzik &
        Markram, 1998.
        The tmgsyn is a dynamic three state implicit resource synapse model.
        The response onset is instantaneous and the decay is exponential.
        It combines a frequency dependent depression and a facilitation
        mechanism that both depend on independent time constants.
        The synaptic targets are chosen by proximity, that is the target pool
        are the cells in closest proximity.

        Parameters
        ----------
        pre_pop - gennetwork.Population
            The presynaptic population
        post_pop - gennetwork.Population
            the postsynaptic population
        target_pool - int
            the number of cells in the target pool
        target_segs - str
            the name of the segments that are possible synaptic targets at the
            postsynaptic population
        divergence - int
            divergence in absolute terms, that is the number of synapses each
            presynaptic cell forms
        tau_1 - numeric
            the time constant of synaptic decay. conforms to the transition
            from the active to the inactive resource state. units of time as in
            neuron standard units
        tau_facil - numeric
            the time constant of facilitation decay. this essentially creates
            the frequency dependence. set to 0 for no facilitation.
        U - numeric
            maximum of postsynaptic response. has to be considered together
            with the weight from the netcon.
        tau_rec - numeric
            time constant of recovery from inactive for recovered state.
            gives depression since inactive resources do not contribute to
            postsynaptic signal. set to 0 for no depression.
        e - numeric
            equilibrium potential of the postsynaptic conductance
        thr - numeric
            threshold for synaptic event at the presynaptic source
        delay - numeric
            delay between presynaptic signal and onset of postsynaptic signal
        weight - numeric
            weight for the netcon object connecting source and target

        Returns
        -------
        None

        Use Cases
        ---------
        >>> tmgsynConnection(nw.population[0], nw.population[1],
                             3, 'prox', 1, 6.0, 0, 0.04, 0, 0, 10, 3, 0)
        A non-facilitating, non-depressing excitatory connection.

        """
        self.init_parameters = locals()
        self.pre_pop = pre_pop
        self.post_pop = post_pop
        pre_pop.add_connection(self)
        post_pop.add_connection(self)
        pre_pop_rad = (np.arange(pre_pop.get_cell_number(), dtype=float) /
                       pre_pop.get_cell_number()) * (2*np.pi)
        post_pop_rad = (np.arange(post_pop.get_cell_number(), dtype=float) /
                        post_pop.get_cell_number()) * (2*np.pi)

        pre_pop_pos = pos(pre_pop_rad)
        post_pop_pos = pos(post_pop_rad)
        pre_cell_target = []
        synapses = []
        netcons = []

        # Setup the Gaussian distribution
        loc = post_pop.get_cell_number() / 2
        gauss = stats.expon(loc=0, scale=scale)
        pdf = gauss.pdf(np.arange(post_pop.get_cell_number()))
        pdf = pdf/pdf.sum()

        for idx, curr_cell_pos in enumerate(pre_pop_pos):

            curr_dist = []
            for post_cell_pos in post_pop_pos:
                curr_dist.append(euclidian_dist(curr_cell_pos, post_cell_pos))

            sort_idc = np.argsort(curr_dist)
            picked_cells = np.random.choice(sort_idc, divergence,
                                            replace=True, p=pdf)
            pre_cell_target.append(picked_cells)
            for target_cell in picked_cells:

                curr_syns = []
                curr_netcons = []

                curr_seg_pool = post_pop[target_cell].get_segs_by_name(target_segs)
                chosen_seg = np.random.choice(curr_seg_pool)
                for seg in chosen_seg:
                    curr_syn = h.tmgsyn(chosen_seg(0.5))
                    curr_syn.tau_1 = tau_1
                    curr_syn.tau_facil = tau_facil
                    curr_syn.U = U
                    curr_syn.e = e
                    curr_syn.tau_rec = tau_rec
                    curr_syns.append(curr_syn)
                    curr_netcon = h.NetCon(pre_pop[idx].soma(0.5)._ref_v,
                                           curr_syn, thr, delay, weight,
                                           sec=pre_pop[idx].soma)
                    curr_netcons.append(curr_netcon)
                    netcons.append(curr_netcons)
                    synapses.append(curr_syns)

        self.netcons = netcons
        self.pre_cell_targets = np.array(pre_cell_target)
        self.synapses = synapses


class tmgsynConnection_old(GenConnection):

    def __init__(self, pre_pop, post_pop,
                 target_pool, target_segs, divergence,
                 tau_1, tau_facil, U, tau_rec, e, thr, delay, weight):
        """Create a connection with tmgsyn as published by Tsodyks, Pawelzik &
        Markram, 1998.
        The tmgsyn is a dynamic three state implicit resource synapse model.
        The response onset is instantaneous and the decay is exponential.
        It combines a frequency dependent depression and a facilitation
        mechanism that both depend on independent time constants.
        The synaptic targets are chosen by proximity, that is the target pool
        are the cells in closest proximity.

        Parameters
        ----------
        pre_pop - gennetwork.Population
            The presynaptic population
        post_pop - gennetwork.Population
            the postsynaptic population
        target_pool - int
            the number of cells in the target pool
        target_segs - str
            the name of the segments that are possible synaptic targets at the
            postsynaptic population
        divergence - int
            divergence in absolute terms, that is the number of synapses each
            presynaptic cell forms
        tau_1 - numeric
            the time constant of synaptic decay. conforms to the transition
            from the active to the inactive resource state. units of time as in
            neuron standard units
        tau_facil - numeric
            the time constant of facilitation decay. this essentially creates
            the frequency dependence. set to 0 for no facilitation.
        U - numeric
            maximum of postsynaptic response. has to be considered together
            with the weight from the netcon.
        tau_rec - numeric
            time constant of recovery from inactive for recovered state.
            gives depression since inactive resources do not contribute to
            postsynaptic signal. set to 0 for no depression.
        e - numeric
            equilibrium potential of the postsynaptic conductance
        thr - numeric
            threshold for synaptic event at the presynaptic source
        delay - numeric
            delay between presynaptic signal and onset of postsynaptic signal
        weight - numeric
            weight for the netcon object connecting source and target

        Returns
        -------
        None

        Use Cases
        ---------
        >>> tmgsynConnection(nw.population[0], nw.population[1],
                             3, 'prox', 1, 6.0, 0, 0.04, 0, 0, 10, 3, 0)
        A non-facilitating, non-depressing excitatory connection.

        """
        self.init_parameters = locals()
        self.pre_pop = pre_pop
        self.post_pop = post_pop
        pre_pop.add_connection(self)
        post_pop.add_connection(self)
        pre_pop_rad = (np.arange(pre_pop.get_cell_number(), dtype=float) /
                       pre_pop.get_cell_number()) * (2*np.pi)
        post_pop_rad = (np.arange(post_pop.get_cell_number(), dtype=float) /
                        post_pop.get_cell_number()) * (2*np.pi)

        pre_pop_pos = pos(pre_pop_rad)
        post_pop_pos = pos(post_pop_rad)
        pre_cell_target = []
        synapses = []
        netcons = []

        for idx, curr_cell_pos in enumerate(pre_pop_pos):
            curr_dist = []
            for post_cell_pos in post_pop_pos:
                curr_dist.append(euclidian_dist(curr_cell_pos, post_cell_pos))

            sort_idc = np.argsort(curr_dist)
            closest_cells = sort_idc[0:target_pool]
            picked_cells = np.random.choice(closest_cells,
                                            divergence,
                                            replace=False)
            pre_cell_target.append(picked_cells)
            for tar_c in picked_cells:

                curr_syns = []
                curr_netcons = []

                curr_seg_pool = post_pop[tar_c].get_segs_by_name(target_segs)
                chosen_seg = np.random.choice(curr_seg_pool)
                for seg in chosen_seg:
                    curr_syn = h.tmgsyn(chosen_seg(0.5))
                    curr_syn.tau_1 = tau_1
                    curr_syn.tau_facil = tau_facil
                    curr_syn.U = U
                    curr_syn.e = e
                    curr_syn.tau_rec = tau_rec
                    curr_syns.append(curr_syn)
                    curr_netcon = h.NetCon(pre_pop[idx].soma(0.5)._ref_v,
                                           curr_syn, thr, delay, weight,
                                           sec=pre_pop[idx].soma)
                    curr_netcons.append(curr_netcon)
                    netcons.append(curr_netcons)
                    synapses.append(curr_syns)

        self.netcons = netcons
        self.pre_cell_targets = np.array(pre_cell_target)
        self.synapses = synapses


class Exp2SynConnection(GenConnection):
    """
    This class connects a pre and a post synaptic population with a Exp2Syn
    synapse.
    """

    def __init__(self, pre_pop, post_pop, target_pool, target_segs, divergence,
                 tau1, tau2, e, thr, delay, weight):
        """
        divergence,
                 tau1, tau2, e, g_max, thr, delay, weight, name = "GC->MC"
                 """
        self.init_parameters = locals()
        self.pre_pop = pre_pop
        self.post_pop = post_pop
        pre_pop.add_connection(self)
        post_pop.add_connection(self)
        pre_pop_rad = (np.arange(pre_pop.get_cell_number(), dtype=float) /
                       pre_pop.get_cell_number()) * (2*np.pi)
        post_pop_rad = (np.arange(post_pop.get_cell_number(), dtype=float) /
                        post_pop.get_cell_number()) * (2*np.pi)
        self.pre_pop_rad = pre_pop_rad
        self.post_pop_rad = post_pop_rad

        pre_pop_pos = pos(pre_pop_rad)
        post_pop_pos = pos(post_pop_rad)
        pre_cell_target = []
        synapses = []
        netcons = []

        for idx, curr_cell_pos in enumerate(pre_pop_pos):

            curr_dist = []
            for post_cell_pos in post_pop_pos:
                curr_dist.append(euclidian_dist(curr_cell_pos, post_cell_pos))

            sort_idc = np.argsort(curr_dist)
            closest_cells = sort_idc[0:target_pool]
            picked_cells = np.random.choice(closest_cells,
                                            divergence,
                                            replace=False)
            pre_cell_target.append(picked_cells)
            for tar_c in picked_cells:
                curr_syns = []
                curr_netcons = []

                curr_seg_pool = post_pop[tar_c].get_segs_by_name(target_segs)
                chosen_seg = np.random.choice(curr_seg_pool)
                for seg in chosen_seg:
                    curr_syn = h.Exp2Syn(chosen_seg(0.5))
                    curr_syn.tau1 = tau1
                    curr_syn.tau2 = tau2
                    curr_syn.e = e
                    curr_syns.append(curr_syn)
                    curr_netcon = h.NetCon(pre_pop[idx].soma(0.5)._ref_v,
                                           curr_syn, thr, delay, weight,
                                           sec=pre_pop[idx].soma)
                    curr_netcons.append(curr_netcon)
                    netcons.append(curr_netcons)
                    synapses.append(curr_syns)

        self.netcons = netcons
        self.pre_cell_targets = np.array(pre_cell_target)
        self.synapses = synapses


class PerforantPathStimulation(object):
    """
    This class connects a pre and a post synaptic population with a Exp2Syn
    synapse.
    """

    def __init__(self, stim, post_pop, n_targets, target_segs,
                 tau1, tau2, e, thr, delay, weight):
        """
        divergence,
                 tau1, tau2, e, g_max, thr, delay, weight, name = "GC->MC"
        """

        self.pre_pop = stim
        self.post_pop = post_pop
        post_pop.add_connection(self)
        synapses = []
        netcons = []

        if type(n_targets) == int:
            # Select n_targets from post_pop
            target_cells = np.random.choice(post_pop.cells, n_targets,
                                            replace=False)
        else:
            target_cells = post_pop.cells[n_targets]

        for curr_cell in target_cells:
            curr_seg_pool = curr_cell.get_segs_by_name(target_segs)
            for seg in curr_seg_pool:
                curr_syn = h.Exp2Syn(seg(0.5))
                curr_syn.tau1 = tau1
                curr_syn.tau2 = tau2
                curr_syn.e = e
                curr_netcon = h.NetCon(stim, curr_syn, thr, delay, weight)
                netcons.append(curr_netcon)
                synapses.append(curr_syn)

        self.netcons = netcons
        self.pre_cell_targets = np.array(target_cells)
        self.synapses = synapses


class PerforantPathPoissonStimulation(object):
    """
    Patterned Perforant Path stimulation as in Yim et al. 2015.
    uses vecevent.mod -> h.VecStim
    """
    def __init__(self, post_pop, t_pattern, spat_pattern, target_segs,
                 tau1, tau2, e, weight):

        post_pop.add_connection(self)
        synapses = []
        netcons = []
        conductances = []

        target_cells = post_pop.cells[spat_pattern]
        self.pre_pop = "Implicit"
        self.vecstim = h.VecStim()
        self.pattern_vec = h.Vector(t_pattern)
        self.vecstim.play(self.pattern_vec)

        for curr_cell in target_cells:
            curr_seg_pool = curr_cell.get_segs_by_name(target_segs)
            curr_conductances = []
            for seg in curr_seg_pool:
                curr_syn = h.Exp2Syn(seg(0.5))
                curr_syn.tau1 = tau1
                curr_syn.tau2 = tau2
                curr_syn.e = e
                curr_netcon = h.NetCon(self.vecstim, curr_syn)
                curr_gvec = h.Vector()
                curr_gvec.record(curr_syn._ref_g)
                curr_conductances.append(curr_gvec)
                curr_netcon.weight[0] = weight
                netcons.append(curr_netcon)
                synapses.append(curr_syn)
                """for event in pattern:
                    curr_netcon.event(event)"""
            conductances.append(curr_conductances)

        self.netcons = netcons
        self.pre_cell_targets = np.array(target_cells)
        self.synapses = synapses
        self.conductances = conductances


class PerforantPathPoissonTmgsyn(GenConnection):
    """
    Patterned Perforant Path simulation as in Yim et al. 2015.
    uses vecevent.mod -> h.VecStim
    """
    def __init__(self, post_pop, t_pattern, spat_pattern, target_segs,
                 tau_1, tau_facil, U, tau_rec, e, weight):

        self.init_parameters = locals()
        post_pop.add_connection(self)
        synapses = []
        netcons = []
        t_pattern = list(t_pattern)  # nrn does not like np.ndarrays?
        target_cells = post_pop[spat_pattern]
        self.pre_pop = 'Implicit'
        self.post_pop = post_pop
        self.vecstim = h.VecStim()
        self.pattern_vec = h.Vector(t_pattern)
        self.vecstim.play(self.pattern_vec)
        conductances = []

        for curr_cell in target_cells:
            curr_seg_pool = curr_cell.get_segs_by_name(target_segs)
            curr_conductances = []
            for seg in curr_seg_pool:
                curr_syn = h.tmgsyn(seg(0.5))
                curr_syn.tau_1 = tau_1
                curr_syn.tau_facil = tau_facil
                curr_syn.U = U
                curr_syn.tau_rec = tau_rec
                curr_syn.e = e
                curr_netcon = h.NetCon(self.vecstim, curr_syn)
                curr_gvec = h.Vector()
                curr_gvec.record(curr_syn._ref_g)
                curr_conductances.append(curr_gvec)
                curr_netcon.weight[0] = weight
                netcons.append(curr_netcon)
                synapses.append(curr_syn)
            conductances.append(curr_conductances)

        self.conductances = conductances
        self.netcons = netcons
        self.pre_cell_targets = np.array(spat_pattern)
        self.synapses = synapses

"""Population ONLY REMAINS IN gennetwork TO KEEP pyDentate RUNNING. THE NEW
IMPLEMENTATION OF POPULATION IS IN genpopulation"""
class Population(object):
    """This is the model of a generic population.
    A population is a number of cells of a specific type derived from
    genneuron.GenNeuron. The GenPopulation object keeps track of all
    incoming and outgoing connections. It is recommended to create Populations
    through the GenNetwork.mk_population interface of a network the population
    is part of.

    Attributes
    ----------
    parent_network - gennetwork.GenNetwork or derived instances
        The network the population takes part in
    cell_type - genneuron.GenNeuron class or subclass thereof
        The cell type making up the population
    cells - list of genneuron.GenNeuron instances
        A list of cells that currently exist within the population
    connections - list of Connection objects
        A list of outgoing and incoming connections

    Methods
    -------
    __init__
    make_cells
    get_cell_number
    record_aps
    plot_aps
    write_aps
    current_clamp_rnd
    current_clamp_range
    voltage_recording
    add_connection

    Use cases
    ---------
    >>> nw = GenNetwork()
    >>> nw.mk_population(GranuleCell, 500)
    Create an empty network and create a population of 500 granule cells in the
    network.
    """

    def __init__(self, cell_type=None, n_cells=None, parent_network=None):
        self.parent_network = parent_network
        self.cell_type = cell_type
        self.cells = []
        self.connections = []
        self.VClamps = []
        self.VClamps_i = []
        self.VRecords = []
        if cell_type and n_cells:
            self.make_cells(cell_type, n_cells)
        self.i = 0

    def SEClamp(self, cells, dur1=200, amp1=0, rs=0.001):
        for x in cells:
            clamp = self.cells[x]._SEClamp(dur1=dur1, amp1=amp1, rs=rs)
            self.VClamps.append(clamp)
            curr_vec = h.Vector()
            curr_vec.record(clamp._ref_i)
            self.VClamps_i.append(curr_vec)

    def voltage_recording(self, cells):
        for x in cells:
            record = self.cells[x]._voltage_recording()
            self.VRecords.append(record)

    def make_cells(self, cell_type, n_cells):
        """Create cells of a certain type

        Parameters
        ----------
        cell_type - genneuron.GenNeuron class of subclass thereof
            the type of the cells to be created
        n_cells - numeric
            number of cells to be created

        Returns
        -------
        None

        Use Cases
        ---------
        >>> popul = Population(parent_network = nw)
        >>> popul.make_cells(GranuleCell, 500)
        Create an empty population within nw and then create 500 granule cells
        """

        if hasattr(self, 'cell_type'):
            if self.cell_type != cell_type:
                raise TypeError("cell_type inconsistent with population")
            else:
                self.cell_type = cell_type

        if not hasattr(self, 'cells'):
            self.cells = []

        for _ in range(n_cells):
            self.cells.append(cell_type())

        self.cells = np.array(self.cells, dtype=object)

    def get_cell_number(self):
        """Return the number of cells"""
        return len(self.cells)

    def record_aps(self):
        counters = [cell._AP_counter() for cell in self.cells]
        self.ap_counters = counters
        return counters

    def plot_aps(self, color='k'):
        cells = []
        for x in self.ap_counters:
            # as_numpy() doesn't work on windows 10 ???
            try:
                cells.append(x[0].as_numpy())
            except:
                cells.append(np.array(x[0]))

        # Workaround for matplotlib bug. plt.eventplot throws error when first
        # element empty
        if not np.array(cells[0]).any():
            cells[0] = np.array([0], dtype=float)

        plt.eventplot(cells, linewidth=2, color=color)

    def write_aps(self, directory='', fname=''):
        if not fname:
            time_tup = time.gmtime()
            time_str = time.asctime(time_tup)
            time_str = '_'.join(time_str.split(' '))
            nw_name = self.parent_network.__class__.name
            pop_name = self.cell_type.name
            fname = nw_name + '_' + pop_name + '_' + time_str
            fname = fname.replace(':', '-')
        if not directory:
            directory = os.getcwd()
        if not os.path.isdir(directory):
            os.mkdir(directory)
        path = directory + '\\' + fname + '.npz'
        try:
            ap_list = [x[0].as_numpy() for x in self.ap_counters]
        except:
            ap_list = [np.array(x[0]) for x in self.ap_counters]
        np.savez(path, *ap_list)

    def perc_active_cells(self):
        try:
            # as_numpy doesn't work on windows 10 ???
            timing_arrays = [x[0].as_numpy() for x in self.ap_counters]
        except:
            timing_arrays = [np.array(x[0]) for x in self.ap_counters]
        active_counter = 0
        for x in timing_arrays:
            if x.size != 0:
                active_counter += 1

        return (active_counter / float(self.get_cell_number())) * 100

    def mk_current_clamp(self, cells, amp=0.3, dur=5, delays=3):
        if not hasattr(cells, '__iter__'):
            cells = np.random.choice(self.get_cell_number(), cells,
                                     replace=False)

        if not hasattr(delays, '__iter__'):
            delays = np.array(delays)

        for cell in cells:
            for delay in delays:
                self.cells[cell]._current_clamp_soma(amp=amp, dur=dur,
                                                     delay=delay)

    def current_clamp_rnd(self, n_cells, amp=0.3, dur=5, delay=3):
        """DEPRECATE"""
        chosen_cells = np.random.choice(self.cells, n_cells, replace=False)

        for x in chosen_cells:
            for y in delay:
                x._current_clamp_soma(amp=amp, dur=dur, delay=y)

        return chosen_cells

    def current_clamp_range(self, n_cells, amp=0.3, dur=5, delay=3):
        """DEPRECATE"""
        if type(n_cells) == int:
            n_cells = range(n_cells)

        for cell in n_cells:
            self.cells[cell]._current_clamp_soma(amp=amp, dur=dur, delay=delay)

    """def voltage_recording(self, cell_type):
        rnd_int = random.randint(0, len(self.cells) - 1)
        soma_v_vec = self.cells[rnd_int]._voltage_recording()
        return soma_v_vec"""

    def add_connection(self, conn):
        self.connections.append(conn)

    def get_properties(self):
        """Get the properties of the network"""
        try:
            ap_time_stamps = [x[0].as_numpy() for x in self.ap_counters]
        except:
            ap_time_stamps = [np.array(x[0]) for x in self.ap_counters]
        ap_numbers = [x[1].n for x in self.ap_counters]
        try:
            v_rec = [x.as_numpy() for x in self.VRecords]
            vclamp_i = [x.as_numpy() for x in self.VClamps_i]
        except:
            v_rec = [np.array(x) for x in self.VRecords]
            vclamp_i = [np.array(x) for x in self.VClamps_i]
        properties = {'parent_network': str(self.parent_network),
                      'cell_type': self.cell_type.name,
                      'cell_number': self.get_cell_number(),
                      'connections': [conn.get_properties()
                                      for conn in self.connections],
                      'ap_time_stamps': ap_time_stamps,
                      'ap_number': ap_numbers,
                      'v_records': v_rec,
                      'VClamps_i': vclamp_i}

        properties
        return properties

    def __str__(self):
        return self.cell_type.name + 'Population'

    def __iter__(self):
        return self

    def __getitem__(self, item):
        return self.cells[item]

    def __next__(self):
        if self.i < (len(self.cells)):
            i = self.i
            self.i += 1
            return self.cells[i]
        else:
            self.i = 0
            raise StopIteration()

    def next(self):
        return self.__next__()

# HELPERS
def pos(rad):
    """
    (x,y) position of a point on a circle with axis origin at (0,0)
    and radius 1.
    x = cx + r * cos(rad) -> x = cos(rad)
    y = cy + r * sin(rad) -> y = sin(rad)

    Returns a list of tuples that give the point of each radian passed.
    """
    x_arr = list(np.cos(rad))
    y_arr = list(np.sin(rad))

    return [(x_arr[idx], y_arr[idx]) for idx in range(len(x_arr))]


def euclidian_dist(p1, p2):
    """ p1 and p2 must both be of len 2 where p1 = (x1,y1); p2 = (x2,y2)"""
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)
