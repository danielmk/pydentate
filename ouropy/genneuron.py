# -*- coding: utf-8 -*-
"""
This module implements the GenNetwork class, which implements generic network
logic.
"""

from neuron import h
from ouropy.gendendrite import GenDendrite
import numpy as np


class GenNeuron(object):
    """This is the model of a generic neuron.

    Attributes
    ----------
    soma - nrn.Section (Default None)
        The soma
    dendrites - list (Default empty)
        A list of gendendrite.GenDendrite
    all_secs - list (Default empty)
        A list of all sections

    Methods
    -------
    __init__
    mk_soma
    mk_dendrite
    conn_post

    Use cases
    ---------
    >>> myNeuron = GenNeuron()
    >>> myNeuron.mk_soma()
    >>> myNeuron.mk_dendrite()
    Ball-and-stick neuron with default geometry
    """

    def mk_soma(self, diam=None, L=None, name=None):
        """Assignes self.soma a hoc section with dimensions diam and L.
        Uses nrn defaults when None. Name defaults to 'soma'.
        Before mk_soma is called, self.soma = None.

        Parameters
        ----------
        diam - numeric (Default from NEURON)
            diameter of the soma
        L - numeric (Default from NEURON)
            length of the soma
        name - str (Default 'soma')
            name of the section

        Returns
        -------
        None

        Use cases
        ---------
        >>> self.mk_soma()
        self.soma becomes section with default values
        >>> self.mk_soma(name = 'myFirstSoma', diam = 5, L = 100)
        self.soma becomes section with name = 'myFirstSoma', diam = 5 and
            L = 100
        """

        if not name:
            name = 'soma'

        if hasattr(self, 'soma'):
            self.all_secs.remove(self.soma)

        self.soma = h.Section(name=name)

        if diam:
            self.soma.diam = diam
        if L:
            self.soma.L = L

        if not hasattr(self, 'all_secs'):
            self.all_secs = []
        self.all_secs.append(self.soma)

    def mk_dendrite(self, n_secs=1, dend_name=None, sec_names=None, diam=None,
                    L=None, soma_loc=1):
        """Adds a dendrite to list self.dendrites. Before first call
        self.dendrites = None. On first call self.dendrite becomes list with 1
        dendrite. Automatically connects the first section of the dendrite to
        the soma. Raises error if self.soma = None.

        Parameters
        ----------
        n_secs - int
            number of sections on the dendrite
        dend_name - str
            the name
        sec_names - list of str
            the names of the sections in order
        diam - list of numerics
            the diameters of the sections
        L - list of numerics
            the diameters of the sections

        Returns
        -------
        None

        Use cases
        ---------
        >>> self.mk_dendrite()
        Create a dendrite with 1 section and default values
        >>> self.mk_dendrite(4, 'myDendrite',
                             ['prox1', 'prox2', 'dist1', 'dist2'],
                             [5,3,3,5], [100,50,50,50])
        Create a dendrite with 4 sections named 'myDendrite' and custom
        geometry and section naming.
        """
        if not hasattr(self, 'dendrites'):
            self.dendrites = []

        if not self.soma:
            raise StandardError("No soma created yet.")

        curr_dend = GenDendrite(dend_name, n_secs, sec_names, diam, L)
        curr_dend.conn_soma(self.soma, soma_loc)
        self.dendrites.append(curr_dend)
        for x in curr_dend:
            self.all_secs.append(x)

    def get_segs_by_name(self, name):
        """Returns a list of sections whose .name matches the name parameter.

        Parameters
        ----------
        name - str or list of str
            The names to gets

        Returns
        -------
        result - list
            List of sections whose .name attribute matches the name parameter

        Use cases
        ---------
        >>> self.get_segs_by_name('proxd')
        Returns segments named 'proxd'
        """

        if 'all' in name:
            return list(self.all_secs)

        result = []
        if type(name) == str:
            for x in self.all_secs:
                if x.name() == name:
                    result.append(x)
        else:
            for x in name:
                if type(x) != str:
                    raise TypeError("All elements of name must be str")
                for y in self.all_secs:
                    if y.name() == x:
                        result.append(y)

        return np.array(result, dtype=np.dtype(object))

    def insert_mechs(self, parameters):
        """Inserts the parameters into the section of the cells.
        See ouropy.parameters

        Parameters
        ----------
        parameters - ouropy.parameters.Parameter or ParameterSet
            A parameter or parameterset contains the mechanism, the segment and
            the value to assign.

        Returns
        -------
        None

        Use Cases
        ---------
        >>> import ouropy.parameters
        >>> myParameters = ouropy.parameters.read_parameters(filename)
        >>> self.insert_mechs(myParameters)
        Insert the parameters loaded from filename into self. See
        ouropy.parameters for details.
        """
        mechanisms = parameters.get_mechs()

        for x in mechanisms.keys():
            sections = self.get_segs_by_name(x)
            for y in sections:
                for z in mechanisms[x]:
                    y.insert(z)

        for x in parameters:
            sections = self.get_segs_by_name(x.sec_name)
            for y in sections:
                setattr(y, x.mech_name, x.value)

    def _current_clamp_soma(self, amp=0.3, dur=500, delay=500):
        """Setup a current clamp at the soma(0.5) section.

        Parameters
        ----------
        amp - numeric
            amplitude of current injection in neurons unit of current (nA)
        dur - numeric
            duration of current injection in neurons unit of time (ms)
        delay - numeric
            start of current injection in neurons unit of time (ms)

        Returns
        -------
        stim - h.IClamp
            the current clamp point process

        Use Cases
        ---------
        >>> self._current_clamp_soma()
        Setup a default current clamp.
        >>> volt_vector, t_vector = self._current_clamp_soma()
        Setup a default current clamp and assign a voltage vector recording at
        soma(0.5) and a time vector.
        """
        if not hasattr(self, 'stim'):
            self.stim = []

        stim = h.IClamp(self.soma(0.5))
        stim.amp = amp     # Too high amps crash the simulation w/o error!
        stim.dur = dur
        stim.delay = delay
        self.stim.append(stim)

        return stim

    def _voltage_recording(self):
        """Record the voltage at the soma of GenNeuron or a subclass.

        Parameters
        ----------
        None

        Returns
        -------
        soma_v_vec - h.Vector
            the current clamp point process

        Use Cases
        ---------
        >>> myNeuron = GenNeuron()
        >>> myNeuron.mk_soma()
        >>> myNeuron._voltage_recording
        """

        soma_v_vec = h.Vector()
        t_vec = h.Vector()
        soma_v_vec.record(self.soma(0.5)._ref_v)
        t_vec.record(h._ref_t)

        return soma_v_vec

    def _AP_counter(self, thr=0):
        """Action Potential counter at the soma of GenNeuron or a subclass.

        Parameters
        ----------
        thr - numeric
            the threshold value for action potential detection

        Returns
        -------
        time_vec - h.Vector
            the current clamp point process
        ap - h.APCount
            the ap counter object that counts number and records time stamps

        Use Cases
        ---------
        >>> myNeuron = GenNeuron()
        >>> myNeuron.mk_soma()
        >>> myNeuron._AP_counter
        """

        ap = h.APCount(self.soma(0.5))
        self.ap = ap
        self.ap.thresh = thr
        self.time_vec = h.Vector()
        self.ap.record(self.time_vec)
        return self.time_vec, self.ap

    def _current_ladder(self, currents, start, stim_dur, stim_interval):
        for idx, x in enumerate(currents):
            delay = (idx * stim_interval) + start
            self._current_clamp_soma(amp=x, dur=stim_dur, delay=delay)

    def _SEClamp(self, dur1=200, amp1=0, rs=0.001):
        self.vclamp = h.SEClamp(self.soma(0.5))
        self.vclamp.dur1 = dur1
        self.vclamp.amp1 = amp1
        self.vclamp.rs = rs
        return self.vclamp

    def __iter__(self):
        return self

    def __next__(self):
        if not hasattr(self, 'i'):
            self.i = 0
        if self.i < (len(self.all_secs)):
            i = self.i
            self.i += 1
            return self.all_secs[i]
        else:
            self.i = 0
            raise StopIteration()

    def next(self):
        return self.__next__()