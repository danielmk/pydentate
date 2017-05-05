# -*- coding: utf-8 -*-
"""
Created on Fri May 05 13:28:44 2017

@author: DanielM
"""

from mpi4py import MPI
from neuron import h

pc = h.ParallelContext()

id = int(pc.id())
nhost = int(pc.nhost())

print "I am", id, "of", nhost