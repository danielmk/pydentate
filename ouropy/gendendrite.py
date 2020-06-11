# -*- coding: utf-8 -*-
"""
This module implements the GenDendrite class, which implements generic dendrite
logic.
"""
from neuron import h


class GenDendrite(object):
    """This is the model of a generic dendrite.

    Attributes:
    name - String (Default None)
        Name of the dendrite
    secs - list (Default None)
        List of the dendrites sections
    soma - nrn.Section (Default None)
        The soma to which the dendrite connects

    Methods:
    __init__
    mk_secs
    conn_soma
    set_diam
    set_L

    Use cases:
    >>> myDend = GenDendrite()
    Creates a default dendrite without sections
    >>> myDend = GenDendrite('myTestDend',4,['prox1','prox2','dist1','dist_2'],
                             [5,5,8,8], [50,50,70,100])
    Creates a dendrite with 4 sections and specified geometry

    """
    def __init__(self, dend_name=None, n_secs=None, sec_names=None, diam=None,
                 L=None):
        self.name = dend_name
        self.secs = None
        self.soma = None
        self._i = 0

        if n_secs:
            self.mk_secs(n_secs, sec_names)

        if diam:
            self.set_diam(diam)

        if L:
            self.set_L(L)

    def mk_secs(self, n_secs=1, sec_names=[]):
        """Makes sections AND connects them. This is because a dendrite is by
        definition made up of connected sections. sec_names has to be a list
        of section names with len = n_secs. If sec_names = None the section
        names are 'sec' + str(number).
        """

        self.secs = []

        if sec_names:
            if not (hasattr(sec_names, '__getitem__')):
                raise TypeError("sec_names should be list or None")
            if len(sec_names) != n_secs:
                raise ValueError("The len of sec_names must equal n_secs")

        for curr_n in range(n_secs):
            if sec_names:
                self.secs.append(h.Section(name=sec_names[curr_n]))
            else:
                self.secs.append(h.Section(name='sec' + str(curr_n)))
            if curr_n > 0:
                self.secs[curr_n].connect(self.secs[curr_n - 1](1))

    def conn_soma(self, soma, soma_loc=1):
        """Connect a soma to the dendrite"""
        if self.soma:
            raise StandardError("Soma already connected")

        if not self.secs:
            raise StandardError("Dendrite has no sections")
        self.soma = soma
        self.secs[0].connect(self.soma(soma_loc))

    def set_diam(self, diam):
        """Change the diameter of the dendrite"""
        if not bool(self.secs):
            raise StandardError("Can't set diameter before sections are made")
        if hasattr(diam, '__getitem__'):
            if len(diam) != len(self.secs):
                raise StandardError("List of diams does not fit n_secs")
            for idx, curr_seg in enumerate(self.secs):
                curr_seg.diam = diam[idx]
            return
        else:
            for curr_seg in self.secs:
                curr_seg.diam = diam

    def set_L(self, L):
        """Change the length of the dendrite"""
        if not bool(self.secs):
            raise Warning("Can't set L before segments are made")
        if hasattr(L, '__getitem__'):
            if len(L) != len(self.secs):
                raise Warning("List of diams does not fit number of segments")
            for idx, curr_seg in enumerate(self.secs):
                curr_seg.L = L[idx]
            return
        else:
            for curr_seg in self.secs:
                curr_seg.L = L

    def __iter__(self):
        return self

    def __next__(self):
        if not self.secs:
            raise StandardError("No sections created yet")
        if self._i < (len(self.secs)):
            i = self._i
            self._i += 1
            return self.secs[i]
        else:
            self._i = 0
            raise StopIteration()

    def next(self):
        return self.__next__()

    def __getitem__(self, key):
        if type(key) == int:
            return self.secs[key]
        else:
            for x in self.secs:
                if x.name == key:
                    return x
            raise KeyError('Key not found')

    def __len__(self):
        return len(self.secs)
